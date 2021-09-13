// #include "multi_filament_header.hpp"
// #include "FilamentJacobianSolver.hpp"
// #include "c_array_functions.hpp"
// #include "FCMfunctions.hpp"
// #include "print_functions.hpp"
// #include "file_functions.hpp"
// #include "filament_initialisation_functions.hpp"
// #include "profilers.hpp"
// #include "CollisionBarrierFilament.hpp"
// #include "spring_link_functions.hpp"
// #include <limits>
// #include <cstddef>

// #include <mpi.h>
// #include <omp.h>
#include <iostream>
// #include <cmath>
// #include <vector>
// #include <fstream>
// #include <stdlib.h>
// #include <stdio.h>


// int Np = 30;// Overwritten later
// #if verbose
//   #include <cassert>     // for DEBUGGING
// #endif

// Simulation parameters are set in "config.hpp"

int main(int argc, char **argv){

using namespace std;

cout << "hello" << endl;

  // cout << "A" << endl;
  //   int totalnodes, myrank;
  //   MPI_Init(&argc, &argv);
  //   MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  //   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //
  //   // Filament setup ==========================================================
  //   // Declare array of Nsw filaments formed of Nworm beads.
  //   // Positions will be kept track of only by rank 0.
  //   vector<Filament> filaments;
  //
  // // Declare linked list struct with all variables needed for collision
  // // barrier linked list. Then initialise.
  // cout << Nsw << " " << max_Nworm << " "<<Nsw*max_Nworm << endl;
  // linked_list linked_list;
  // cout << myrank << " " << totalnodes << endl;
  // linked_list.initialise(myrank, totalnodes, filaments);
  // cout << "X";

	return 0;
}
