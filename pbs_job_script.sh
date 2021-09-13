#PBS -lwalltime=0:01:00
#PBS -lselect=1:ncpus=24:mem=120gb
#PBS -o /rds/general/user/akt12/home/MultiFilamentFCM-Adam/pbs_output/
#PBS -e /rds/general/user/akt12/home/MultiFilamentFCM-Adam/pbs_output/
continue_from_file = xxCONTINUEFROM # leave blank if desired
module load fix_unwritable_tmp
cd $PBS_O_WORKDIR/ # cd to directory script was called from
source run_wrapper_cx2.sh # Load libraries
mpiexec ./MPIFilament $continue_from_file |& tee ./pbs_output/$PBS_JOBID.txt # http://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/computing/high-throughput-computing/configuring-mpi-jobs/
