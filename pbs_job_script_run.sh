#PBS -lwalltime=23:59:00
#PBS -lselect=6:ncpus=24:mem=32gb
#PBS -o /Users/adam/Documents/GitHub/force-coupling-method-codes/MultiFilamentFCM-Adam/pbs_output/
#PBS -e /Users/adam/Documents/GitHub/force-coupling-method-codes/MultiFilamentFCM-Adam/pbs_output/
module load fix_unwritable_tmp
cd $PBS_O_WORKDIR/ # cd to directory script was called from
source run_wrapper_cx2.sh # Load libraries
cd 1912102117/Test-Nnetfil1-Kap3p6-NPTSx144y144z144/
mpiexec ./MPIFilament Test-Nnetfil1-Kap3p6-NPTSx144y144z144  |& tee ../../pbs_output/$PBS_JOBID.Test-Nnetfil1-Kap3p6-NPTSx144y144z144.txt # http://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/computing/high-throughput-computing/configuring-mpi-jobs/
cd ../..
