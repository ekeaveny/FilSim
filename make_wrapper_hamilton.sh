#!/bin/bash
module purge
module load intel/2019.5 intelmpi/intel/2019.6
module load fftw/impi/intel/2.1.5

module load intel-suite mpi
module load armadillo

make MPI_FIL_hamilton
echo "--------------------------------------------------"
echo "    MPIFilament compiled and libraries loaded."
source run_wrapper_hamilton.sh
