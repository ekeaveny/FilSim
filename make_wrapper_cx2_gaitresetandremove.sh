#!/bin/bash
module load intel-suite mpi
module load armadillo
module load fftw/2.1.5-double
make MPI_FIL_cx2_gaitresetandremove
echo "--------------------------------------------------"
echo "    MPIFilament compiled and libraries loaded."
source run_wrapper_cx2.sh
