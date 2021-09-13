#PBS -lwalltime=xxWALLTIME
#PBS -lselect=xxNUMNODES:ncpus=xxNUMCPUS:mem=xxMEMORY:mpiprocs=xxMPIPROCS
#PBS -o xxPWD/pbs_output/
#PBS -e xxPWD/pbs_output/
module load fix_unwritable_tmp
cd $PBS_O_WORKDIR/ # cd to directory script was called from
source run_wrapper_cx2.sh # Load libraries
cd xxTIMESTAMP/xxSIM/
mpiexec ./MPIFilament xxCONTINUEFROM xxRESTARTFROMTIMEZERO |& tee ../../pbs_output/$PBS_JOBID.xxSIM.txt # http://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/computing/high-throughput-computing/configuring-mpi-jobs/
cd ../..
