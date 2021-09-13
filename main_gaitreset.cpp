#include "multi_filament_header.hpp"
#include "FilamentJacobianSolver.hpp"
#include "c_array_functions.hpp"
#include "FCMfunctions.hpp"
#include "print_functions.hpp"
#include "file_functions.hpp"
#include "filament_initialisation_functions.hpp"
#include "profilers.hpp"
#include "CollisionBarrierFilament.hpp"
#include "spring_link_functions.hpp"
#include <limits>
#include <cstddef>

int Np = Nsw*Nworm; // Overwritten later
#if verbose
  #include <cassert>     // for DEBUGGING
#endif

// Simulation parameters are set in "config.hpp"

int main(int argc, char **argv){

	int totalnodes, myrank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


	// For FFT code to work best, NPTS_X, NPTS_Y and NPTS_Z should be
	// 2^N*P for some low prime (or low prime power) p. It also works best if
	// NPTS_X is divisible by the number of nodes, as this is the direction
	// that FFTW splits the domain for parallel processing. None of these things
	// are necessary, however.                               -- [AKT, 27/08/19]
	/*
	   if (NPTS_X%totalnodes!=0 || NPTS_Y%totalnodes!=0 || NPTS_Z%totalnodes!=0) {
	    printf("Parallel processing/FCM error: For FCM code to work, "
	           "NPTS_X, NPTS_Y and NPTS_Z must be divisible by the total number of "
	           "nodes requested. "
	           "You have requested NPTS_X = %ld, NPTS_Y = %ld, NPTS_Z = %ld, totalnodes = %i. "
	           "\n",
	           NPTS_X, NPTS_Y, NPTS_Z, totalnodes);
	    assert(NPTS_X%totalnodes==0);
	    assert(NPTS_Y%totalnodes==0);
	    assert(NPTS_Z%totalnodes==0);
	   }
	 */

	// Filament setup ==========================================================
	// Declare empty array of filaments
	// Positions will be kept track of only by rank 0.
	// Individual filaments must be initialised and pushed onto this array
	// using (e.g.)
	//     "Filament new_filament(Nworm);"
	//     "filaments.push_back(new_filament);"
	// within filament initialisation functions.
	// (This is a change from mid-2020 so older initialisation functions
	// might segfault if they haven't been changed to include this.)
	vector<Filament> filaments;

	//cout << Nsw << " " << Nworm << endl;
	//filaments[0].print_properties();

	// Not used. KAP is multiplied by this
	vec BendingFactors(Nsw,fill::ones); //

	// BeadNumbers is used in saving the parameter file later.
	// Filament lengths are set to Nworm by default.
	ivec BeadNumbers(Nsw,fill::ones);
	BeadNumbers *= Nworm;

	#if EnableSpringLinks
		// Declare network spring links
		spring_links spring_links;
	#endif

	bool checkpointing = false;
	bool restart_using_checkpoint_as_initialisation = false;
	if(argc==2) {
		// Continue from checkpoint as if nothing has happened
		// (i.e. write to old files)
		// argv[1] is checkpoint file
		checkpointing = true;
	} else if (argc==3) {
		// Start from nt=0 using checkpoint .bak file as initialisation only
		// (i.e. write to new files)
		// argv[1] is checkpoint file
		// argv[2] is anything, e.g. 'restart', to toggle this mode.
		restart_using_checkpoint_as_initialisation = true;
	}

	// Output setup information
	if(myrank == 0) {
		// Note Np here will be the value before it's overwritten later.
		print_parameter_values_in_table(Np, LfcmBox_x, LfcmBox_y, LfcmBox_z,
		                                L, omega, dt);
		set_precision_of_screen_output(15);
	}

	#if PROFILING
		profiler profiler;
	#endif

	// Position initialisation (read in from checkpoint if asked)
	int nt_start = 1; // Starting timestep
	int forceArrayLengthFilled = 0; // } so that filament knows where it
	int stateArrayLengthFilled = 0; // } sits in the force/state vectors
	if(myrank == 0) {
		if(checkpointing) {
			// Checkpointing. Read README.md about checkpointing limitations
			// Load nt, X, q, U, Xt, qt, Ut, lam, lam1, lam2 from file and
			// into  filaments .
			// NOTE Can not cope with filaments with different swim phases.
			tie(nt_start, forceArrayLengthFilled, stateArrayLengthFilled)
			        = continue_from_checkpoint(argv[1],filaments);
			// Apply filament properties (e.g. bending modulus) as if a new
			// simulation. This should be fine.
			apply_filament_properties(filaments);

			cout << "Filament initialisation from checkpoint ("
			     << argv[1] << ".bak) complete. " << endl;

			// Spring links between the ends of the network filaments
			#if EnableSpringLinks
				spring_links.read_spring_links_from_file(argv[1]);
				cout << spring_links.size() << " spring links read in "
				     << "from data file (" << argv[1] << "-springlinks.dat)."
				     << endl;
			#endif

			cout <<"Continuing from timestep " << nt_start << "/" << TimeSteps
			     << "." <<endl << endl;
		}
		else if (restart_using_checkpoint_as_initialisation) {
			// Read in as in checkpointing, but set nt_start = 0.
			tie(nt_start, forceArrayLengthFilled, stateArrayLengthFilled)
			        = continue_from_checkpoint(argv[1],filaments);
			nt_start = 1; // Overwrite

			apply_filament_properties(filaments);

			cout << "Filament initialisation from checkpoint ("
			     << argv[1] << ".bak) complete. " << endl;

			// Spring links between the ends of the network filaments
			#if EnableSpringLinks
				spring_links.read_spring_links_from_file(argv[1]);
				cout << spring_links.size() << " spring links read in "
				     << "from data file (" << argv[1] << "-springlinks.dat)."
				     << endl;
			#endif

			cout << "Starting from timestep 1/" << TimeSteps
			     << "." <<endl << endl;

			// Save parameter file
			for (int i=0; i<filaments.size(); i++) {
				BeadNumbers[i] = filaments[i].length();
			}
			save_parameter_values_to_file(SimulationConfigName,
			                              LfcmBox_x, LfcmBox_y, LfcmBox_z,
			                              BeadNumbers, BendingFactors,
			                              filaments);

			// Write initial configuration to file
			save_data_column_names_to_file(SimulationDataName, filaments);
			save_data_to_file(SimulationDataName,filaments,0);
			save_backup_data_to_file(SimulationBackupDataName,filaments,0);
			if (SaveExtendedSwimmerData) {
				save_data_column_names_to_swimmer_velocity_file(SimulationSwimmerDataName, filaments);
				save_swimmer_velocity_data_to_file(SimulationSwimmerDataName,filaments,NULL,NULL,NULL,0);
                save_data_column_names_to_swimmer_collision_file(SimulationSwimmerDataCollisionName, filaments);
				save_swimmer_collision_data_to_file(SimulationSwimmerDataCollisionName,filaments,NULL,0);
			}
			#if GaitReset
				save_data_column_names_to_swimmer_velocity_file(SimulationGaitResetName, filaments);
				save_swimmer_velocity_data_to_file(SimulationGaitResetName,filaments,NULL,NULL,NULL,0);
			#endif
		}
		else {
			// We are not using a checkpoint
			cout << "No checkpoint used; filament initialisation from fresh." << endl;
			// Initialise positions. Note this edits  filaments  in-place.
			#if EnableSpringLinks
				cout << "Initialising filament positions with spring links..." << endl;
				tie(forceArrayLengthFilled, stateArrayLengthFilled)
				        = filament_position_initialisation(filaments,
				                                           spring_links);
			#else
				cout << "Initialising filament positions..." << endl;
				tie(forceArrayLengthFilled, stateArrayLengthFilled)
				        = filament_position_initialisation(filaments);
			#endif
			cout << "Initialising filament properties..." << endl;
			apply_filament_properties(filaments);

			// Start a new data file.
			save_data_column_names_to_file(SimulationDataName, filaments);

			cout<<"Filament initialisation complete."<<endl<<endl;

			// Save parameter file
			for (int i=0; i<filaments.size(); i++) {
				BeadNumbers[i] = filaments[i].length();
			}
			save_parameter_values_to_file(SimulationConfigName,
			                              LfcmBox_x, LfcmBox_y, LfcmBox_z,
			                              BeadNumbers, BendingFactors,
			                              filaments);

			// Write initial configuration to file
			save_data_to_file(SimulationDataName,filaments,0);
			save_backup_data_to_file(SimulationBackupDataName,filaments,0);
			if (SaveExtendedSwimmerData) {
				save_data_column_names_to_swimmer_velocity_file(SimulationSwimmerDataName, filaments);
				save_swimmer_velocity_data_to_file(SimulationSwimmerDataName,filaments,NULL,NULL,NULL,0);
                save_data_column_names_to_swimmer_collision_file(SimulationSwimmerDataCollisionName, filaments);
                save_swimmer_collision_data_to_file(SimulationSwimmerDataCollisionName,filaments,NULL,0);
			}
			#if GaitReset
				save_data_column_names_to_swimmer_velocity_file(SimulationGaitResetName, filaments);
				save_swimmer_velocity_data_to_file(SimulationGaitResetName,filaments,NULL,NULL,NULL,0);
			#endif
			// Spring links between the ends of the network filaments
			#if EnableSpringLinks
				spring_links.decide_spring_links(filaments);
				cout << spring_links.size() << " spring links formed." << endl;
				spring_links.save_spring_links_to_file(SimulationSpringLinkDataName);
			#endif
		}

		#if verbose
			filaments[0].printEverything();
			std::cin.get();
		#endif

		// #if PROFILING
		// 	profiler.totalTime.start();
		// #endif
	}

	// compute global number of particles on node 0, then broadcast to all.
	Np = 0;
	for(int nn = 0; nn<filaments.size(); ++nn) {
		Np += filaments[nn].length();
	}
	MPI_Bcast(&Np, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Initialise custom solver
	filamentJsolver Jsolver(Np);
	bool check;
	int current_plot_step = 0;

	// dynamic allocation (heap) for large matrices.
	vec Error(6*Np,fill::zeros);
	vec NewError(6*Np,fill::zeros);
	vec Update(6*Np,fill::zeros);

	// Broadcast nt_start
	MPI_Bcast(&nt_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Declare linked list struct with all variables needed for collision
	// barrier linked list. Then initialise.
	LinkedList linked_list;
	linked_list.initialise(myrank, totalnodes, filaments);

	// Declare fcm struct with all variables needed for the force coupling
	// method. Then initialise.
	Fcm fcm;

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point 000A]" << endl;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	fcm.initialise(myrank);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point 00A]" << endl;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	// First swimmer's centre of mass
	// Used later to measure speed of the interesting filament
	int interesting_filament;
	vec com0_0, com_0(3);
	vec com0_1, com_1(3);
	vec com0_2, com_2(3);
	if(myrank == 0) {
		interesting_filament = 0;
		com0_0 = filaments[interesting_filament].getCOM(); // CoM at time t=0
		com_0 = com0_0;// Current CoM
		com0_1 = filaments[interesting_filament].getCOM(); // CoM at time t=0
		com_1 = com0_1;// Current CoM
		com0_2 = filaments[interesting_filament].getCOM(); // CoM at time t=0
		com_2 = com0_1;// Current CoM
	}

	// Fixed width of timestep/Broyden's iter counter
	int len_TimeSteps = to_string(TimeSteps).length();
	int len_broyden_maxiter = to_string(broyden_maxiter).length();

	int iter; // Broyden's iteration counter
	string temp_str;
	// Time stepping begins ====================================================
	for (int nt=nt_start; nt<=timesteps; ++nt) {
		//cout << "[Point 0A]" << myrank << endl;
		//MPI_Barrier(MPI_COMM_WORLD); // remove later

		// Run twice for gait reset.
		// gait_reset_i = 2: normal step
		// gait_reset_i = 1: fake step; no collisions for the swimmer
		// gait_reset_i = 0: fake step; no collisions at all
		for (int gait_reset_i=0; gait_reset_i<=2; ++gait_reset_i)
		{

			if(myrank == 0) {
				// Print column headings
				if ((nt-nt_start) % 20 == 0) {
					cout << "[ " << SimulationName << " ]" << endl;
					temp_str = "Timestep";
					temp_str.resize(len_TimeSteps*2+1,' ');
					cout << " " << temp_str << "  ";
					temp_str = "Broyden its";
					temp_str.resize(len_broyden_maxiter,' ');
					cout << "(" << temp_str << ") ";
					#if PROFILING
						profiler.print_names();
					#endif
					cout << endl;
				}

				// Initially, we update forces and torques within each filament.
				// For the collision barrier and fluid solve, we copy them into
				// fcm.F, fcm.T,
				// such that we can easily plug in different fluid solvers
				// and collision barrier functions that are agnostic of swimmers.

				//cout << "[Point A]" << myrank << endl;
				#if PROFILING
					profiler.totalTimeFCM.reset();
					profiler.totalTimeCopy.reset();
					profiler.totalTimeBuildJacobian.reset();
					profiler.totalTimeErrorCheck.reset();
					profiler.totalTimeSolveJacobian.reset();
					profiler.totalTimestepTime.reset();
					profiler.totalTimeCollisionBarrier.reset();

					profiler.fcm_make_zero.reset();
					profiler.fcm_gaussian_setup.reset();
					profiler.fcm_force_distribution.reset();
					profiler.fcm_fft_forward.reset();
					profiler.fcm_fft_backward.reset();
					profiler.fcm_flow_solve.reset();
					profiler.fcm_particle_velocities_rotations.reset();

					profiler.totalTimestepTime.start();
					profiler.totalTime.start();

					profiler.totalTimeCopy.start();
				#endif

				// #pragma omp parallel for // Probably not worth parallelising
				for(int nn=0; nn < filaments.size(); ++nn) {
					//cout << nn << " / " << filaments.size() << endl;
					filaments[nn].initialGuess(); // Xs = 2*X - Xt etc...
					filaments[nn].unRobotArmify();
					filaments[nn].setZeroForcesTorques();
					// filaments[nn].setZeroLambdas(); // Gait reset special. Delete after.
					// if(nn==0){
					// 	filaments[nn].printEverything();
					// 	cin.get();
					// }
					// Driving (for swimming) is handled by applyElasticTorques.
					// Make sure SwimmingHelixAlpha and SwimmingHelixBeta are set.
					filaments[nn].applyElasticTorques(nt);
					//filaments[nn].printEverything();
					//cin.get();
					// if(nn==0){
					// 	filaments[nn].printEverything();
					// 	cin.get();
					// }
					filaments[nn].applyConstraintForcesTorques();
					//filaments[nn].printEverything();
					//cin.get();
					// if(nn==0){
					// 	filaments[nn].printEverything();
					// 	cin.get();
					// }
					filaments[nn].applyExternalForcesTorques(nt);
					//filaments[nn].printEverything();
					//cin.get();
					// if(nn==0){
					// 	filaments[nn].printEverything();
					// 	cin.get();
					// }
					// cout << "END OF OUTPUTS" << endl;
				}

				#if EnableSpringLinks
					// Apply spring link forces
					spring_links.apply_spring_link_forces(filaments);
				#endif

				// Collision barrier force added later
				//cout << "[Point B]" << myrank << endl;
				fcm.assign_filament_data(filaments);

				#if PROFILING
					profiler.totalTimeCopy.end();
					profiler.totalTimeBuildJacobian.start();
				#endif
				//cout << "[Point C]" << myrank << endl;
				// analytical approximate Jacobian for each filament.
				Jsolver.buildJacobian(filaments,nt);

				#if PROFILING
					profiler.totalTimeBuildJacobian.end();
				#endif

				#if MPIverbose
					fcm.print_data("After forces/torques applied [Start]",
					               "After forces/torques applied [End]",true);
				#endif
			}
			//cout << "[Point Cb]" << myrank << endl;
			//MPI_Barrier(MPI_COMM_WORLD); // remove later
			//cout << "[Point D]" << myrank << endl;
			// MPI communicate FCM.Y, .F, .T to everyone
			fcm.mpi_broadcast();
			//MPI_Barrier(MPI_COMM_WORLD); // remove later
			//cout << "[Point E]" << myrank << endl;
			//MPI_Barrier(MPI_COMM_WORLD); // remove later

			if (gait_reset_i>=1) {

				if (myrank == 0) {
					#if PROFILING
						profiler.totalTimeCollisionBarrier.start();
					#endif

					// Set up linked list for given filament positions
					linked_list.link(filaments);
					//cout << "[Point F]" << myrank << endl;
				}
				//MPI_Barrier(MPI_COMM_WORLD); // remove later
				//cout << "[Point G]" << myrank << endl;

				linked_list.mpi_broadcast();
				//MPI_Barrier(MPI_COMM_WORLD); // remove later
				//cout << "[Point H]" << myrank << endl;

                if (gait_reset_i == 2){
                    linked_list.apply_collision_barrier(fcm,nt);
                } else {
                    linked_list.apply_collision_barrier_not_swimmer(fcm,nt);
                }
				//MPI_Barrier(MPI_COMM_WORLD); // remove later
				//cout << "[Point I]" << myrank << endl;

				if (myrank == 0) {
					#if PROFILING
						profiler.totalTimeCollisionBarrier.end();
					#endif
				}
				//MPI_Barrier(MPI_COMM_WORLD); // remove later
				//cout << "[Point J]" << myrank << endl;
			}
			/*
			   MPI_Barrier(MPI_COMM_WORLD);
			   if (myrank == 0) {
			    for (int i = 0; i < Np; i++) {
			        for (int j = 0; j < 3; j++) {
			            cout << "F["<<i<<"]["<<j<<"] = " << fcm.F[i][j] << endl;
			        }
			    }
			    for (int i = 0; i < Np; i++) {
			        for (int j = 0; j < 3; j++) {
			            cout << "T["<<i<<"]["<<j<<"] = " << fcm.T[i][j] << endl;
			        }
			    }
			   }
			   cin.get();
			 */

			if(nt>initialFrictionSteps) {
				#if PROFILING
					if(myrank == 0) {
						profiler.totalTimeFCM.start();
					}
				#endif
				//MPI_Barrier(MPI_COMM_WORLD); // remove later
				//cout << "[Point K]" << myrank << endl;

				#if PROFILING
					fcm.mobility_solve(&profiler);
				#else
					fcm.mobility_solve();
				#endif
				//MPI_Barrier(MPI_COMM_WORLD); // remove later
				//cout << "[Point L]" << myrank << endl;
				#if PROFILING
					if (myrank == 0) {
						profiler.totalTimeFCM.end();
					}
				#endif
			} else {
				//cout << "[Point KL2]" << myrank << endl;
				fcm.friction_mobility_solve();
			}

			/*
			   MPI_Barrier(MPI_COMM_WORLD);
			   if (myrank == 0) {
			    for (int i = 0; i < Np; i++) {
			        for (int j = 0; j < 3; j++) {
			            cout << "V["<<i<<"]["<<j<<"] = " << fcm.V[i][j] << endl;
			        }
			    }
			   }
			   cin.get();
			 */

			if(myrank == 0) {
				//cout << "[Point M]" << myrank << endl;
				#if MPIverbose
					fcm.print_data("After velocities found [Start]",
					               "After velocities found [End]",false);
				#endif
				#if PROFILING
					profiler.totalTimeErrorCheck.start();
				#endif
				//cout << "[Point N]" << myrank << endl;
				// integration and constraint error
				check = ErrorCheckFCM(filaments, Error, fcm.V, fcm.W, nt);
				//cout << "[Point O]" << myrank << endl;
				#if flush_on
					cout << "\r" << flush;
				#else
					cout << " " << setw(len_TimeSteps) << nt << "/"
					     << TimeSteps << "  ";
				#endif
				#if PROFILING
					profiler.totalTimeErrorCheck.end();
				#endif

				iter = 0;
				// double DampingAlpha = 1.;
			}

			MPI_Bcast(&iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
			//cout << "[Point P]" << myrank << endl;
			MPI_Bcast(&check, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
			//cout << "[Point Q]" << myrank << endl;

			while (check && iter<broyden_maxiter) {
				if(myrank == 0) {
					#if flush_on
						cout << " " << setw(len_TimeSteps) << nt << "/"
						     << TimeSteps << "  ";
						cout << "(" << setw(len_broyden_maxiter) << iter << ") ";
					#endif
					#if PROFILING
						profiler.totalTimeSolveJacobian.start();
					#endif
					// update Jacobian of solver (would be needed for good Broyden)
					// JacobianLU.assignCurrentStep();

					// solve J*Update = Error
					//cout << "[Point R]" << myrank << endl;
					Jsolver.solve(filaments, Error, Update, iter);
					//cout << "[Point S]" << myrank << endl;
					// cout << "a MPI J solve. my rank = " << myrank << endl;
					#if PROFILING
						profiler.totalTimeSolveJacobian.end();
						profiler.totalTimeCopy.start();
					#endif
					// damping
					// Update *= DampingAlpha;

					#pragma omp parallel for
					for(int nn=0; nn<filaments.size(); ++nn) {
						filaments[nn].update(Update); // apply update to each fil.
						filaments[nn].unRobotArmify();
						filaments[nn].setZeroForcesTorques();
						// Driving (for swimming) is handled by applyElasticTorques.
						// Make sure SwimmingHelixAlpha & SwimmingHelixBeta are set.
						filaments[nn].applyElasticTorques(nt);
						filaments[nn].applyConstraintForcesTorques();
						filaments[nn].applyExternalForcesTorques(nt);
					}

					#if EnableSpringLinks
						// Apply spring link forces
						spring_links.apply_spring_link_forces(filaments);
					#endif

					// Collision barrier force added later
					//cout << "[Point T]" << myrank << endl;
					fcm.assign_filament_data(filaments);

					#if PROFILING
						profiler.totalTimeCopy.end();
					#endif
				}
				/// MPI communicate Y, F, T to everyone
				//cout << "[Point U]" << myrank << endl;
				fcm.mpi_broadcast();

				if(gait_reset_i >= 1) {

					if (myrank == 0) {
						//cout << "[Point V]" << myrank << endl;
						#if PROFILING
							profiler.totalTimeCollisionBarrier.start();
						#endif
					}

					//cout << "[Point W]" << myrank << endl;
                    if(gait_reset_i == 2){
                        linked_list.apply_collision_barrier(fcm,nt);
                    } else {
                        linked_list.apply_collision_barrier_not_swimmer(fcm,nt);
                    }


					/*  if (myrank == 0) {
					    cout << "B" << endl;
					    for (int i=0; i<Np; ++i) {
					        for (int j=0; j<3; ++j) {
					            //fcm.F[i][j] += linked_list.F_total[i][j];
					            //cout << "F_total["<<i<<"]["<<j<<"] = " << linked_list.F_total[i][j] << " ==> F["<<i<<"]["<<j<<"] = " << fcm.F[i][j] << endl;
					            cout << "F["<<i<<"]["<<j<<"] = " << fcm.F[i][j] << endl;
					        }
					    }
					   }
					   //cin.get();
					 */


					if (myrank == 0) {
						#if PROFILING
							profiler.totalTimeCollisionBarrier.end();
						#endif
					}
				}

				// fluid solve
				if(nt>initialFrictionSteps) {
					#if PROFILING
						profiler.totalTimeFCM.start();
					#endif
					//cout << "[Point X]" << myrank << endl;
					#if PROFILING
						fcm.mobility_solve(&profiler);
					#else
						fcm.mobility_solve();
					#endif
					#if PROFILING
						profiler.totalTimeFCM.end();
					#endif
				} else {
					//cout << "[Point X2]" << myrank << endl;
					fcm.friction_mobility_solve();
				}


				if(myrank == 0) {
					#if MPIverbose
						fcm.print_data("Before error check [Start]",
						               "Before error check [End]",false);
					#endif

					#if PROFILING
						profiler.totalTimeErrorCheck.start();
					#endif
					//cout << "[Point Y]" << myrank << endl;
					check = ErrorCheckFCM(filaments, NewError, fcm.V, fcm.W, nt);
					//cout << "[Point Z]" << myrank << endl;
					#if PROFILING
						profiler.totalTimeErrorCheck.end();
					#endif

					// Damping
					// if(NewError.norm() > Error.norm()){
					//   DampingAlpha /= 2.;
					// }
					// else{
					//   DampingAlpha = .5;
					// }

					#if verbose
						cout << "ERROR" <<endl<<Error <<endl;
						cout << "state after error check\n";
						filaments[0].printEverything(); std::cin.get();
						cout << "===========" << endl;
						std::cin.get();
					#endif

					// update for Bad Broyden's
					vec tmp(Nbroy,fill::zeros);
					vec DeltaError(Nbroy);
					DeltaError = NewError - Error;
					double DeltaErrorNorm = norm(DeltaError);

					#if PROFILING
						profiler.totalTimeSolveJacobian.start();
					#endif
					//cout << "[Point ZA]" << myrank << endl;
					Jsolver.solve(filaments, DeltaError, tmp, iter);
					#if PROFILING
						profiler.totalTimeSolveJacobian.end();
					#endif
					// cout << "b MPI J solve. my rank = " << myrank << endl;

					//cout << "[Point ZB]" << myrank << endl;
					tmp *= -1./DeltaErrorNorm;
					tmp -= Update/DeltaErrorNorm;
					Jsolver.addCmatCol(tmp, iter);
					tmp = DeltaError/DeltaErrorNorm;
					Jsolver.addDmatCol(tmp, iter);

					Error = NewError;

					#if verbose
						cout << endl << endl << "nt " << nt << ", iter " << iter
						     <<" end of Newton step \n" << endl;
					#endif

					#if flush_on
						cout << ""<< '\r' << std::flush;
					#endif
					++iter;

					#if verbose
						if(iter==broyden_maxiter) {
							cout << "Newton iteration not converged \n";
							std::cin.get();
						}
					#endif

				}
				//cout << "[Point ZC]" << myrank << endl;
				MPI_Bcast(&iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&check, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
				//cout << "[Point ZD]" << myrank << endl;
			}

			// MPI_Barrier(MPI_COMM_WORLD); // remove later
			if(myrank == 0) {
				#if verbose
					cout << endl << "exited Newton iterations \n \n";
					cout << "Error\n"<< Error << endl;
				#endif
				#if MPIverbose
					fcm.print_data("After step completed [Start]",
					               "After step completed [End]",true);
				#endif
				#if verbose
					filaments[0].printEverything(); std::cin.get();
				#endif

				// End of timestep printout
				vec average_v;
				if (gait_reset_i==0) {
					com_0 = filaments[interesting_filament].getCOM();
					average_v = (com_0-com0_0)/(dt);
					com0_0 = com_0;
				} else if (gait_reset_i == 1) {
					com_1 = filaments[interesting_filament].getCOM();
					average_v = (com_1-com0_1)/(dt);
					com0_1 = com_1;
				} else {
                    com_2 = filaments[interesting_filament].getCOM();
                    average_v = (com_2-com0_2)/(dt);
                    com0_2 = com_1;
                }

				if(gait_reset_i == 2) {

					// shift state variables back by one time step
					#pragma omp parallel for
					for(int nn=0; nn<filaments.size(); ++nn) {
						filaments[nn].step(nt);
					}

					#if flush_on
						cout << " " << setw(len_TimeSteps) << nt << "/"
						     << TimeSteps << "  ";
						cout << "(" << setw(len_broyden_maxiter) << iter << ") ";
					#else
						cout << "(" << setw(len_broyden_maxiter) << iter << ") ";
					#endif

					++current_plot_step;
					if (current_plot_step==plot_steps) {
						save_data_to_file(SimulationDataName,filaments,nt);
						save_backup_data_to_file(SimulationBackupDataName,filaments,nt);
						current_plot_step = 0;
					}
					if (SaveExtendedSwimmerData) { // Saves every timestep. You may want to change this.
						save_swimmer_velocity_data_to_file(SimulationSwimmerDataName,
						                                   filaments,
						                                   fcm.V, fcm.W, fcm.F,
						                                   nt);
						save_swimmer_collision_data_to_file(SimulationSwimmerDataCollisionName,
						                                    filaments,
						                                    linked_list.F_total,
						                                    nt);
					}
				} else {
					cout << "(" << setw(len_broyden_maxiter) << iter << ") ";

					save_swimmer_velocity_data_to_file(SimulationGaitResetName,
					                                   filaments,
					                                   fcm.V, fcm.W, fcm.F,
					                                   nt);
				}

				#if PROFILING
					profiler.totalTimestepTime.end();
					profiler.totalTime.end();
					profiler.print_all(nt);
				#endif

				// THIS WORKS, I'VE JUST COMMENTED IT OUT FOR NEATER OUTPUT
				cout << "[Fil" << interesting_filament << " avgV="
				     << std::setprecision(2)
				     << average_v(0) << "," << average_v(1) << "," << average_v(2)
				     << "]";
				// if(linked_list.min_dist2 < 6.25) {//4.41){
				// 	cout << "[Min dist="
				// 	     << std::setprecision(3) << std::fixed << sqrt(linked_list.min_dist2) << std::scientific
				// 	     << "]";
				// }
				cout << endl;
			}

			MPI_Barrier(MPI_COMM_WORLD); // this one has to stay!

		} // Gait reset end

	}

	if (myrank == 0) {
		#if PROFILING
			profiler.totalTime.print_total();
		#endif
	}

	// Free everything
	fcm.free_memory();

	MPI_Finalize();
	return 0;
}
