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

	vector<Filament> filaments;

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

	// randSeed now set in config

	for (int sep_i = 0; sep_i < 80; sep_i++)
	{
        double sep = 2.0 + 0.1*sep_i;
        int num_trials = 5;
		cout << "===== SEP = " << sep << " =====" << endl;
		for (int trial_i = 0; trial_i < num_trials; trial_i++)
		{
            arma_rng::set_seed(randSeed);

            MPI_Barrier(MPI_COMM_WORLD);
            cout << "[Point P1] " << myrank << endl;

			cout << "---- TRIAL = " << trial_i << " ----" << endl;
			// Position initialisation (read in from checkpoint if asked)
			int nt_start = 1; // Starting timestep
			int forceArrayLengthFilled = 0; // } so that filament knows where it
			int stateArrayLengthFilled = 0; // } sits in the force/state vectors
			int forceArrayLengthFilled2 = 0; // } so that filament knows where it
			int stateArrayLengthFilled2 = 0; // } sits in the force/state vectors

            MPI_Barrier(MPI_COMM_WORLD);
            cout << "[Point P2] " << myrank << endl;

			if(myrank == 0) {
				if(checkpointing) {
                    cout << "Continuing from checkpoint..." << endl;
					// Checkpointing. Read README.md about checkpointing limitations
					// Load nt, X, q, U, Xt, qt, Ut, lam, lam1, lam2 from file and
					// into  filaments .
					// NOTE Can not cope with filaments with different swim phases.
					tie(nt_start, forceArrayLengthFilled, stateArrayLengthFilled)
					        = continue_from_checkpoint(argv[1],filaments);

					cout << "Adding two beads..." << endl;
                    arma_rng::set_seed(randSeed+1+trial_i+sep_i*num_trials);

					tie(forceArrayLengthFilled2, stateArrayLengthFilled2)
					        = two_bead_initialisation(filaments, sep);

					cout << "Two beads added." << endl;

					forceArrayLengthFilled += forceArrayLengthFilled2;
					stateArrayLengthFilled += stateArrayLengthFilled2;

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
                    cout << "Initialising from checkpoint..." << endl;
					// Read in as in checkpointing, but set nt_start = 0.
					tie(nt_start, forceArrayLengthFilled, stateArrayLengthFilled)
					        = continue_from_checkpoint(argv[1],filaments);
					nt_start = 1; // Overwrite

					cout << "Adding two beads..." << endl;
                    arma_rng::set_seed(randSeed+1+trial_i+sep_i*num_trials);

					tie(forceArrayLengthFilled2, stateArrayLengthFilled2)
					        = two_bead_initialisation(filaments, sep);

					cout << "Two beads added." << endl;

					forceArrayLengthFilled += forceArrayLengthFilled2;
					stateArrayLengthFilled += stateArrayLengthFilled2;

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
					save_parameter_values_to_file(SimulationTwoParticleConfigName,
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
					}
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

					cout << "Adding two beads..." << endl;
                    arma_rng::set_seed(randSeed+1+trial_i+sep_i*num_trials);

					tie(forceArrayLengthFilled2, stateArrayLengthFilled2)
					        = two_bead_initialisation(filaments, sep);

					cout << "Two beads added." << endl;

					forceArrayLengthFilled += forceArrayLengthFilled2;
					stateArrayLengthFilled += stateArrayLengthFilled2;

					cout << "Initialising filament properties..." << endl;
					apply_filament_properties(filaments);

					// Start a new data file.
					save_data_column_names_to_file(SimulationDataName, filaments);

					cout<<"Filament initialisation complete."<<endl<<endl;

					// Save parameter file
					for (int i=0; i<filaments.size(); i++) {
						BeadNumbers[i] = filaments[i].length();
					}

					save_parameter_values_to_file(SimulationTwoParticleConfigName,
					                              LfcmBox_x, LfcmBox_y, LfcmBox_z,
					                              BeadNumbers, BendingFactors,
					                              filaments);


					// Write initial configuration to file
					save_data_to_file(SimulationDataName,filaments,0);
					save_backup_data_to_file(SimulationBackupDataName,filaments,0);

					if (SaveExtendedSwimmerData) {
						save_data_column_names_to_swimmer_velocity_file(SimulationSwimmerDataName, filaments);
						save_swimmer_velocity_data_to_file(SimulationSwimmerDataName,filaments,NULL,NULL,NULL,0);
					}

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
				//      profiler.totalTime.start();
				// #endif
			}

            MPI_Barrier(MPI_COMM_WORLD);

			// compute global number of particles on node 0, then broadcast to all.
			Np = 0;
			for(int nn = 0; nn<filaments.size(); ++nn) {
				Np += filaments[nn].length();
			}
			cout << "Np is " << Np << endl;
			MPI_Bcast(&Np, 1, MPI_INT, 0, MPI_COMM_WORLD);

			// Initialise custom solver
			// filamentJsolver Jsolver(Np);
			bool check;
			int current_plot_step = 0;

			// dynamic allocation (heap) for large matrices.
			// vec Error(6*Np,fill::zeros);
			// vec NewError(6*Np,fill::zeros);
			// vec Update(6*Np,fill::zeros);

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
			// int interesting_filament;
			// vec com0, com(3);
			// if(myrank == 0) {
			//      interesting_filament = 0;
			//      com0 = filaments[interesting_filament].getCOM(); // CoM at time t=0
			//      com = com0;// Current CoM
			// }

			// Fixed width of timestep/Broyden's iter counter
			int len_TimeSteps = to_string(TimeSteps).length();
			int len_broyden_maxiter = to_string(broyden_maxiter).length();

			int iter; // Broyden's iteration counter
			string temp_str;
			// Time stepping begins ====================================================
			// for (int nt=nt_start; nt<=timesteps; ++nt) {
			int nt=nt_start;
			//cout << "[Point 0A]" << myrank << endl;
			//MPI_Barrier(MPI_COMM_WORLD); // remove later
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
					//filaments[nn].printEverything();
					//cin.get();
					// Driving (for swimming) is handled by applyElasticTorques.
					// Make sure SwimmingHelixAlpha and SwimmingHelixBeta are set.
					filaments[nn].applyElasticTorques(nt);
					//filaments[nn].printEverything();
					//cin.get();
					filaments[nn].applyConstraintForcesTorques();
					//filaments[nn].printEverything();
					//cin.get();
					filaments[nn].applyExternalForcesTorques(nt);
					//filaments[nn].printEverything();
					//cin.get();
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
				// Jsolver.buildJacobian(filaments,nt);

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

			linked_list.apply_collision_barrier_not_test_particles(fcm,nt);
			//MPI_Barrier(MPI_COMM_WORLD); // remove later
			//cout << "[Point I]" << myrank << endl;

			if (myrank == 0) {
				#if PROFILING
					profiler.totalTimeCollisionBarrier.end();
				#endif
			}
			//MPI_Barrier(MPI_COMM_WORLD); // remove later
			//cout << "[Point J]" << myrank << endl;

			#if PROFILING
				if(myrank == 0) {
					profiler.totalTimeFCM.start();
				}
			#endif
			//MPI_Barrier(MPI_COMM_WORLD); // remove later
			//cout << "[Point K]" << myrank << endl;

			// cout << "N_sw is " << N_sw << endl;

			// cout << "Particle 0, according to f.b.Xs, at "
			//      << filaments[N_sw].beads[0].Xs(0) << ", "
			//      << filaments[N_sw].beads[0].Xs(1) << ", "
			//      << filaments[N_sw].beads[0].Xs(2) << endl;
			// cout << "Particle 1, according to f.b.Xs, at "
			//      << filaments[N_sw+1].beads[0].Xs(0) << ", "
			//      << filaments[N_sw+1].beads[0].Xs(1) << ", "
			//      << filaments[N_sw+1].beads[0].Xs(2) << endl;


			// cout << "Particle 0, according to fcm.Y, at "
			//      << fcm.Y[Np-2][0] << ", "
			//      << fcm.Y[Np-2][1] << ", "
			//      << fcm.Y[Np-2][2] << endl;
			// cout << "Particle 1, according to fcm.Y, at "
			//      << fcm.Y[Np-1][0] << ", "
			//      << fcm.Y[Np-1][1] << ", "
			//      << fcm.Y[Np-1][2] << endl;

			// Two particle test. Place a force on the first bead (the one at zero)

			double t1 = fcm.Y[Np-2][0] - fcm.Y[Np-1][0];
			double t2 = fcm.Y[Np-2][1] - fcm.Y[Np-1][1];
			double t3 = fcm.Y[Np-2][2] - fcm.Y[Np-1][2];

			t1 = t1 - LfcmBox_y * floor(t1/LfcmBox_y + 0.5);
			t2 = t2 - LfcmBox_z * floor(t2/LfcmBox_z + 0.5);
			t3 = t3 - LfcmBox_x * floor(t3/LfcmBox_x + 0.5);

			double tt = sqrt(t1*t1 + t2*t2 + t3*t3);
			t1 = t1/tt;
			t2 = t2/tt;
			t3 = t3/tt;

			// 90 rotation ACW
			double n1 = -t2;
			double n2 = t1;
			double n3 = 0;
			double nn = sqrt(n1*n1 + n2*n2);
			n1 = n1/nn;
			n2 = n2/nn;

			// b = t x n
			double b1 = t2*n3-t3*n2;
			double b2 = t3*n1-t1*n3;
			double b3 = t1*n2-t2*n1;
			double bb = sqrt(b1*b1 + b2*b2 + b3*b3);
			b1 = b1/bb;
			b2 = b2/bb;
			b3 = b3/bb;

			double f = 1.0;
			for(int ii=0; ii<6; ii++) {
				MPI_Barrier(MPI_COMM_WORLD);

                // for(int ll=0; ll<Np; ll++){
                //     for(int mm=0; mm<Np; mm++){
                //         fcm.F[ll][0] = 0;
                //         fcm.T[ll][0] = 0;
                //         fcm.F[ll][1] = 0;
                //         fcm.T[ll][1] = 0;
                //         fcm.F[ll][2] = 0;
                //         fcm.T[ll][2] = 0;
                //     }
                // }

				fcm.F[Np-2][0] = 0;
				fcm.F[Np-2][1] = 0;
				fcm.F[Np-2][2] = 0;
				fcm.F[Np-1][0] = 0;
				fcm.F[Np-1][1] = 0;
				fcm.F[Np-1][2] = 0;
				fcm.T[Np-2][0] = 0;
				fcm.T[Np-2][1] = 0;
				fcm.T[Np-2][2] = 0;
				fcm.T[Np-1][0] = 0;
				fcm.T[Np-1][1] = 0;
				fcm.T[Np-1][2] = 0;

				if(ii == 0) {
					fcm.F[Np-2][0] = f*t1;
					fcm.F[Np-2][1] = f*t2;
					fcm.F[Np-2][2] = f*t3;
				} else if(ii == 1) {
					fcm.F[Np-2][0] = f*n1;
					fcm.F[Np-2][1] = f*n2;
					fcm.F[Np-2][2] = f*n3;
				} else if(ii == 2) {
					fcm.F[Np-2][0] = f*b1;
					fcm.F[Np-2][1] = f*b2;
					fcm.F[Np-2][2] = f*b3;
				} else if(ii == 3) {
					fcm.T[Np-2][0] = f*t1;
					fcm.T[Np-2][1] = f*t2;
					fcm.T[Np-2][2] = f*t3;
				} else if(ii == 4) {
					fcm.T[Np-2][0] = f*n1;
					fcm.T[Np-2][1] = f*n2;
					fcm.T[Np-2][2] = f*n3;
				} else {
					fcm.T[Np-2][0] = f*b1;
					fcm.T[Np-2][1] = f*b2;
					fcm.T[Np-2][2] = f*b3;
				}

				// Add f(=1) to F and T
				// if(ii == 0) {
				//      fcm.F[Np-2][0] += f*t1;
				//      fcm.F[Np-2][1] += f*t2;
				//      fcm.F[Np-2][2] += f*t3;
				// } else if(ii == 1) {
				//      fcm.F[Np-2][0] -= f*t1;
				//      fcm.F[Np-2][1] -= f*t2;
				//      fcm.F[Np-2][2] -= f*t3;
				//
				//      fcm.F[Np-2][0] += f*n1;
				//      fcm.F[Np-2][1] += f*n2;
				//      fcm.F[Np-2][2] += f*n3;
				// } else if(ii == 2) {
				//      fcm.F[Np-2][0] -= f*n1;
				//      fcm.F[Np-2][1] -= f*n2;
				//      fcm.F[Np-2][2] -= f*n3;
				//
				//      fcm.F[Np-2][0] += f*b1;
				//      fcm.F[Np-2][1] += f*b2;
				//      fcm.F[Np-2][2] += f*b3;
				// } else if(ii == 3) {
				//      fcm.F[Np-2][0] -= f*b1;
				//      fcm.F[Np-2][1] -= f*b2;
				//      fcm.F[Np-2][2] -= f*b3;
				//
				//      fcm.T[Np-2][0] += f*t1;
				//      fcm.T[Np-2][1] += f*t2;
				//      fcm.T[Np-2][2] += f*t3;
				// } else if(ii == 4) {
				//      fcm.T[Np-2][0] -= f*t1;
				//      fcm.T[Np-2][1] -= f*t2;
				//      fcm.T[Np-2][2] -= f*t3;
				//
				//      fcm.T[Np-2][0] += f*n1;
				//      fcm.T[Np-2][1] += f*n2;
				//      fcm.T[Np-2][2] += f*n3;
				// } else {
				//      fcm.T[Np-2][0] -= f*n1;
				//      fcm.T[Np-2][1] -= f*n2;
				//      fcm.T[Np-2][2] -= f*n3;
				//
				//      fcm.T[Np-2][0] += f*b1;
				//      fcm.T[Np-2][1] += f*b2;
				//      fcm.T[Np-2][2] += f*b3;
				// }

				MPI_Barrier(MPI_COMM_WORLD);

				#if PROFILING
					fcm.mobility_solve(&profiler);
				#else
					fcm.mobility_solve();
				#endif

				MPI_Barrier(MPI_COMM_WORLD);
				if(myrank == 0) {
					cout << ii << ": done." << endl;
					ofstream OutputFile (SimulationTwoParticleName,ios::app);

					OutputFile.precision(10);
					OutputFile.setf(ios::scientific);
					OutputFile.setf(ios::showpoint);

					OutputFile << sep << " "
					           << trial_i << " "
					           << ii << " "

					           << tt << " " << t1 << " " << t2 << " " << t3 << " "
					        //
					           << fcm.F[Np-2][0]*t1 + fcm.F[Np-2][1]*t2 + fcm.F[Np-2][2]*t3<< " "
					           << fcm.F[Np-2][0]*n1 + fcm.F[Np-2][1]*n2 + fcm.F[Np-2][2]*n3<< " "
					           << fcm.F[Np-2][0]*b1 + fcm.F[Np-2][1]*b2 + fcm.F[Np-2][2]*b3<< " "
					           << fcm.T[Np-2][0]*t1 + fcm.T[Np-2][1]*t2 + fcm.T[Np-2][2]*t3<< " "
					           << fcm.T[Np-2][0]*n1 + fcm.T[Np-2][1]*n2 + fcm.T[Np-2][2]*n3<< " "
					           << fcm.T[Np-2][0]*b1 + fcm.T[Np-2][1]*b2 + fcm.T[Np-2][2]*b3<< " "

					        // << fcm.F[Np-1][0]*t1 + fcm.F[Np-1][1]*t2 + fcm.F[Np-1][2]*t3<< " "
					        // << fcm.F[Np-1][0]*n1 + fcm.F[Np-1][1]*n2 + fcm.F[Np-1][2]*n3<< " "
					        // << fcm.F[Np-1][0]*b1 + fcm.F[Np-1][1]*b2 + fcm.F[Np-1][2]*b3<< " "
					        // << fcm.T[Np-1][0]*t1 + fcm.T[Np-1][1]*t2 + fcm.T[Np-1][2]*t3<< " "
					        // << fcm.T[Np-1][0]*n1 + fcm.T[Np-1][1]*n2 + fcm.T[Np-1][2]*n3<< " "
					        // << fcm.T[Np-1][0]*b1 + fcm.T[Np-1][1]*b2 + fcm.T[Np-1][2]*b3<< " "

					           << fcm.V[Np-2][0]*t1 + fcm.V[Np-2][1]*t2 + fcm.V[Np-2][2]*t3<< " "
					           << fcm.V[Np-2][0]*n1 + fcm.V[Np-2][1]*n2 + fcm.V[Np-2][2]*n3<< " "
					           << fcm.V[Np-2][0]*b1 + fcm.V[Np-2][1]*b2 + fcm.V[Np-2][2]*b3<< " "
					           << fcm.W[Np-2][0]*t1 + fcm.W[Np-2][1]*t2 + fcm.W[Np-2][2]*t3<< " "
					           << fcm.W[Np-2][0]*n1 + fcm.W[Np-2][1]*n2 + fcm.W[Np-2][2]*n3<< " "
					           << fcm.W[Np-2][0]*b1 + fcm.W[Np-2][1]*b2 + fcm.W[Np-2][2]*b3<< " "
					           << fcm.V[Np-1][0]*t1 + fcm.V[Np-1][1]*t2 + fcm.V[Np-1][2]*t3<< " "
					           << fcm.V[Np-1][0]*n1 + fcm.V[Np-1][1]*n2 + fcm.V[Np-1][2]*n3<< " "
					           << fcm.V[Np-1][0]*b1 + fcm.V[Np-1][1]*b2 + fcm.V[Np-1][2]*b3<< " "
					           << fcm.W[Np-1][0]*t1 + fcm.W[Np-1][1]*t2 + fcm.W[Np-1][2]*t3<< " "
					           << fcm.W[Np-1][0]*n1 + fcm.W[Np-1][1]*n2 + fcm.W[Np-1][2]*n3<< " "
					           << fcm.W[Np-1][0]*b1 + fcm.W[Np-1][1]*b2 + fcm.W[Np-1][2]*b3<< " "

					           << fcm.Y[Np-2][0] << " "
					           << fcm.Y[Np-2][1] << " "
					           << fcm.Y[Np-2][2] << " "
					           << fcm.Y[Np-1][0] << " "
					           << fcm.Y[Np-1][1] << " "
					           << fcm.Y[Np-1][2] << " "

                               << fcm.F[Np-2][0] << " "
                               << fcm.F[Np-2][1] << " "
                               << fcm.F[Np-2][2] << " "
                               << fcm.V[Np-2][0] << " "
                               << fcm.V[Np-2][1] << " "
                               << fcm.V[Np-2][2] << " ";

                               for(int pp=0; pp<Np; pp++){
                                   OutputFile
                                   << fcm.F[pp][0] << " "
                                   << fcm.F[pp][1] << " "
                                   << fcm.F[pp][2] << " "
                                   << fcm.V[pp][0] << " "
                                   << fcm.V[pp][1] << " "
                                   << fcm.V[pp][2] << " ";
                               }

					           OutputFile << endl;
					OutputFile.close();
				}
			}
			//MPI_Barrier(MPI_COMM_WORLD); // remove later
			//cout << "[Point L]" << myrank << endl;
			#if PROFILING
				if (myrank == 0) {
					profiler.totalTimeFCM.end();
				}
			#endif


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
				// check = ErrorCheckFCM(filaments, Error, fcm.V, fcm.W, nt);
				//cout << "[Point O]" << myrank << endl;
				#if flush_on
					cout << "\r" << flush;
				#else
					cout << " " << setw(len_TimeSteps) << 0 << "/"
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

			// MPI_Barrier(MPI_COMM_WORLD); // remove later
			if(myrank == 0) {

				#if PROFILING
					profiler.totalTimestepTime.end();
					profiler.totalTime.end();
					profiler.print_all(nt);
				#endif

				cout << endl;
			}

			MPI_Barrier(MPI_COMM_WORLD); // this one has to stay!

// }

			if (myrank == 0) {
				#if PROFILING
					profiler.totalTime.print_total();
				#endif
			}

			MPI_Barrier(MPI_COMM_WORLD); // this one has to stay!

			// Free everything
			fcm.free_memory();

			// Remove last two filaments (the test particles)
			// filaments.pop_back();
			// filaments.pop_back();

			// Empty filaments
			filaments.clear();

		} // End trial loop
	} // End sep loop

	MPI_Barrier(MPI_COMM_WORLD); // this one has to stay!
	MPI_Finalize();

	// cout << "[Point T]" << myrank << endl;
	//
	// cin.get();
	//
	// cout << "[Point U]" << myrank << endl;

	return 0;
}
