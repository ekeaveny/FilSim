#include "multi_filament_header.hpp"
#include "c_array_functions.hpp"
#include "FCMfunctions.hpp"
#include "print_functions.hpp"
#include "file_functions.hpp"
#include "filament_initialisation_functions.hpp"
#include "profilers.hpp"
#include "CollisionBarrierFilament.hpp"
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

	// For FCM code to work, the number of grid points has to be divisible by
	// the total number of nodes requested.
	if (NPTS_X%totalnodes!=0 || NPTS_Y%totalnodes!=0 || NPTS_Z%totalnodes!=0) {
		printf("Parallel processing/FCM error: For FCM code to work, "
		       "NPTS_X, NPTS_Y and NPTS_Z must be divisible by the total number of "
		       "nodes requested. "
		       "You have requested NPTS_X = %i, NPTS_Y = %i, NPTS_Z = %i, totalnodes = %i. "
		       "\n",
		       NPTS_X, NPTS_Y, NPTS_Z, totalnodes);
		assert(NPTS_X%totalnodes==0);
		assert(NPTS_Y%totalnodes==0);
		assert(NPTS_Z%totalnodes==0);
	}

	MPI_Barrier(MPI_COMM_WORLD); // remove later

	// Declare variables for everyone, even though only rank 0 needs them.
	// That's because of scoping.

	// Initialise custom solver
	bool check;
	int current_plot_step = 0;

	// Filament setup ==========================================================
	// Declare array of Nsw filaments. Positions will be kept track of only
	// by rank 0.
	std::vector<Filament> filaments(Nsw);

	// NOTE I don't think this works
	// Set different filament lengths, bending moduli. The bending moduli can
	// be overwritten later.
	vec BendingFactors(Nsw,fill::ones); // KAP are multiplied by these
	vec BeadNumbers(Nsw,fill::ones); // } Filament lengths are set to
	BeadNumbers *= Nworm; // } Nworm by default
	for(int nn = 0; nn < Nsw; ++nn) {
		filaments[nn].set_Kap_N_ah(KAP*BendingFactors[nn],
		                           int(BeadNumbers[nn]), 1.);
	}

	// compute global number of particles
	Np = 0;
	for(int nn = 0; nn<Nsw; ++nn) {
		Np += filaments[nn].length();
	}

	// Output setup information
	if(myrank == 0) {
		print_parameter_values_in_table(Np, LfcmBox_x, LfcmBox_y, LfcmBox_z,
		                                L, omega, dt);
		set_precision_of_screen_output(15);
		save_parameter_values_to_file(SimulationConfigName,
		                              LfcmBox_x, LfcmBox_y, LfcmBox_z,
		                              BeadNumbers, BendingFactors);
	}

	// dynamic allocation (heap) for large matrices.
	vec Error(6*Np,fill::zeros);
	vec NewError(6*Np,fill::zeros);
	vec Update(6*Np,fill::zeros);

	#if PROFILING
		profiler profiler;
	#endif

	// Position initialisation (read in from checkpoint if asked)
	int nt_start = 1; // Starting timestep
	int forceArrayLengthFilled = 0; // } so that filament knows where it
	int stateArrayLengthFilled = 0; // } sits in the force/state vectors
	if(myrank == 0) {
		if(argc>=2) {
			// Checkpointing. Read README.md about checkpointing limitations
			// Load nt, X, q, U, Xt, qt, Ut, lam, lam1, lam2 from file and
			// into  filaments .
			// NOTE Can not cope with filaments of different lengths
			// NOTE Can not cope with filaments with different swim phases.
			tie(nt_start, forceArrayLengthFilled, stateArrayLengthFilled)
			        = continue_from_checkpoint(argv[1],filaments);
			// Apply filament properties (e.g. bending modulus) as if a new
			// simulation. This should be fine.
			apply_filament_properties(filaments);

			cout<<"Filament initialisation from checkpoint ("
			    << argv[1] << ".bak) complete. " << endl
			    <<"Restarting from timestep " << nt_start << "/" << TimeSteps
			    << "." <<endl << endl;
		}
		else {
			// We are not using a checkpoint, so start a new data file.
			save_data_column_names_to_file(SimulationDataName, filaments);
			// Initialise positions. Note this edits  filaments  in-place.
			tie(forceArrayLengthFilled, stateArrayLengthFilled)
			        = filament_position_initialisation(filaments);
			apply_filament_properties(filaments);

			cout<<"Filament initialisation complete."<<endl<<endl;
			// Write initial configuration to file
			save_data_to_file(SimulationDataName,filaments,0);
			save_backup_data_to_file(SimulationBackupDataName,filaments,0);
		}

		#if verbose
			filaments[0].printEverything();
			std::cin.get();
		#endif

		#if PROFILING
			profiler.totalTime.start();
		#endif
	}

	// Declare linked list struct with all variables needed for collision
	// barrier linked list. Then initialise.
	linked_list linked_list;
	linked_list.initialise(myrank, totalnodes, filaments);

	// Declare fcm struct with all variables needed for the force coupling
	// method. Then initialise.
	fcm fcm;
	fcm.initialise(myrank);

	// for (int ii=0; ii<fcm.local_nx; ++ii) {
	// 	for (int jj=0; jj<NPTS_XY; ++jj) {
	// 		for (int kk=0; kk<fcm.pad; ++kk) {
	// 			fcm.fx[ii][jj][kk] = std::sin(ii)*std::sin(jj)*std::sin(kk);
	// 			fcm.fy[ii][jj][kk] = std::sin(ii)*std::sin(jj)*std::sin(kk);
	// 			fcm.fz[ii][jj][kk] = std::sin(ii)*std::sin(jj)*std::sin(kk);
	// 		}
	// 	}
	// }


	// MPI communicate FCM.Y, .F, .T to everyone
	//fcm.mpi_broadcast();


	fcm.mobility_solve();

    // for (int ii=0; ii<fcm.local_nx; ++ii) {
    //     for (int jj=0; jj<NPTS_XY; ++jj) {
    //         for (int kk=0; kk<fcm.pad; ++kk) {
    //             cout << ii << "," << jj << "," << kk << ": "
    //             << fcm.ux[ii][jj][kk] << ","
    //             << fcm.uy[ii][jj][kk] << ","
    //             << fcm.uz[ii][jj][kk] << endl;
    //         }
    //     }
    // }


	// Free everything
	fcm.free_memory();

	MPI_Finalize();
	return 0;
}
