#!/bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/rds/general/user/akt12/home/OpenBLAS
export OPENBLAS_NUM_THREADS=1
module load intel-suite mpi
module load armadillo
module load fftw/2.1.5-double
echo "--------------------------------------------------"
echo "          MPIFilament ready to run."
echo ""
echo " NB: OPENBLAS_NUM_THREADS is set to 1 by default."
echo "     Override this by typing (e.g.)"
echo "     'export OPENBLAS_NUM_THREADS=16'."
echo "--------------------------------------------------"
