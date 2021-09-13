// inclusion guards
#ifndef FILCOLBAR_IS_INCLUDED
#define FILCOLBAR_IS_INCLUDED

//#include "multi_filament_header.hpp"

using namespace arma;
/*
vec VecPeriodic(const vec& XX, const vec& YY);

// Collision barrier repulsion force between beads. Works for periodic box and
// does not include repulsion between neighbouring beads in filament.
// Note: strength F^S must be greater than largest K_B/L^2.
void applyFilamentCollisionBarrierIntrinsic(std::vector<Filament>& filaments);
/*
// Eric's collision barrier function
//void collbarrier(double **Y, double **F, int *map, int *head, int *list, int loc_cell_start, int loc_cell_stop, double Rrefsq, double Fref, double rad);
*/
struct LinkedList {
	int *map, *list, *head;
	int M_X, M_Y, M_Z, ncell, loc_cell_start, loc_cell_stop, mapsize;
	double Rref, Fref;
	double min_dist2 = 999999;

	// Have to initialise with constant lengths.
	// These are declared "on the heap" rather than "on the stack" because otherwise large arrays cause segfaults. (Thanks, C++)
	int(*global_segment_number)[max_Nworm] = new int[Nsw][max_Nworm]; // This is how you declare 2D arrays on the heap
	int* global_segment_to_local_filament_number = new int[Nsw*max_Nworm];
	int* global_segment_to_local_segment_number = new int[Nsw*max_Nworm]; // Change to Nsw and it crashes

	double F_on_node[Nsw*max_Nworm][3] = {{ 0 }};
	double F_total[Nsw*max_Nworm][3] = {{ 0 }};

	void initialise(int myrank, int totalnodes, std::vector<Filament>& filaments);

	int icell(int M_X, int M_Y, int M_Z, int x, int y, int z);

	void bulkmap(int *map, int M_X, int M_Y, int M_Z);

	void apply_collision_barrier(Fcm& fcm, int nt); /*,
	                                  int *map, int *head, int *list,
	                                  int loc_cell_start, int loc_cell_stop,
	                                  int *global_segment_to_local_filament_number,
	                                  int *global_segment_to_local_segment_number); */

	void apply_collision_barrier_not_swimmer(Fcm& fcm, int nt);
	void apply_collision_barrier_not_test_particles(Fcm& fcm, int nt);

	void link(vector<Filament>& filaments);

	void set_global_segment_number(vector<Filament>& filaments); /*,
	                               int global_segment_number[Nsw][max_Nworm],
	                               int global_segment_to_local_filament_number[Nsw*max_Nworm],
	                               int global_segment_to_local_segment_number[Nsw*max_Nworm]); */

	void mpi_broadcast();
};

#endif
