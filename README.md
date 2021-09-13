# README #
AKT 1/6/19

## Compiling/running

### Compile/run a parallel batch on cx2 using qsub
* Compiling:
    1. Set the config data you want to run the batch for in `batch_compile.py`. This includes whether the script should continue from a previous checkpoint. This script uses the config template `config_template.hpp` to create temporary `config.hpp` files, which are used when compiling each item in the batch, and then deleted.
    2. Run `python batch_compile.py`.
    3. Each of the compiled executables and their config files are stored in the folder `[timestamp]/SimName0`. Any required checkpoint files will also be here.

* Running:
    1. **Beware of `screen`**: you might get an `qsub: could not create/open tmp file /var/tmp/pbsscrptwbKVMy for script` error if you are inside a `screen` window.
    2. Set the wall time, number of nodes and number of CPUs/node you want to give the batch in `batch_run.py`. (Note that the total number of CPUs must divide into `NPTS_X`, `NPTS_Y` and `NPTS_Z` in `config.hpp` [check].) This script uses the config template `pbs_job_script_run_template.hpp` to create temporary `pbs_job_script_run.hpp` files, which are used when compiling each item in the batch. We could delete them afterwards if we want, but we don't.
    3. Run `python batch_run.py [timestamp]`, using the timestamp from the compiling process.  
    4. Check stats: `qstat`.
    5. To view screen output as it runs: look in `[timestamp]/SimName0/pbs_output` folder.

### Compile/run a single parallel simulation on cx2 using qsub
* Compiling:
    1. Set the config data you want to run the simulation for in `config.hpp`.
    2. Run `source make_wrapper_cx2.sh`

* Running:
    1. **Beware of `screen`**: you might get an `qsub: could not create/open tmp file /var/tmp/pbsscrptwbKVMy for script` error if you are inside a `screen` window.
    2. Set the wall time, number of nodes and number of CPUs/node you want to give the batch in `pbs_job_script.sh`. (Note that the total number of CPUs must divide into `NPTS_XY` and `NPTS_Z` in `config.hpp`.)
    3. Run `qsub pbs_job_script.sh`.
    4. Check stats: `qstat`.
    5. To view screen output as it runs: look in `pbs_output` folder.    

### Compile/run a single non-parallel simulation on cx2 locally (i.e. not on the cluster)
* Compiling:
    1. Set the config data you want to run the simulation for in `config.hpp`.
    2. Run `source make_wrapper_cx2.sh`

* Running:
    1. If `source make_wrapper_cx2.sh` hasn't been called this session, run `source run_wrapper_cx2.sh`.
    2. Run `./MPIFilament`

### Compile/run a single non-parallel/parallel simulation on gehrig/rizzuto
* Compiling:
    1. `make MPI_FIL`

* Running:
    1. `export OPENBLAS_NUM_THREADS=1`. If you don't do this it will be SLOW. (See below for optimal number)
    2. **FOR NON-PARALLEL**: `./MPIFilament`
    3. **FOR PARALLEL**: `mpiexec -np 4 ./MPIFilament`




## Editing ##

* edit `main.cpp`
    * Worms with different lengths: `BeadNumbers`
    * Worms with different KAP: `BendingFactors`
* config file: `config.hpp`
    * Number of worms
    * (RPY only) Number of nodes: also in `config.hpp` under `myOpenMPthreads`. For FCM, this is controlled automatically by fftw.
    * Periodic box size: change `NPTS_XY` and `NPTS_Z` in `config.hpp` because `LfcmBox_xy = NPTS_XY*dx` and `LfcmBox_z = NPTS_Z*dx` (in `multi_filament_header.hpp`);


## Notes ##

* Swimming? Set `SwimmingHelixBeta=1`
* `OPENBLAS_NUM_THREADS`, set in the `export` statement, is the number of threads used by Armadillo. For short filaments, more than 1 is unlikely to help. For larger filaments, where the filament matrix solve is harder, maybe more will help.
* `myOpenMPthreads` is commented out of `config.hpp`, as it is the number of parallel loops for doing RPY. For FCM, it is irrelevant.
* `omp parallel for` on `main.cpp:265` (ish) slows things down on Gehrig, maybe because 30-odd filaments isn't enough to make it worthwhile. Bear in mind!
* For FFT code to work best, `NPTS_X`, `NPTS_Y` and `NPTS_Z` should be 2^N*P for some low prime (or low prime power) P. This is just how FFTW works. It also works best if `NPTS_X` is divisible by the number of nodes, as this is the direction that FFTW splits the domain for parallel processing. None of these things are necessary, however.     

### Checkpointing ###
* Checkpointing is implemented but naively. It assumes that whatever's in `config.hpp` is also right for the restarted file. It reads in `nt, X, q, U, Xt, qt, Ut, lam, lam1, lam2`.
* Every `plot_steps` timesteps, when the program saves to file, it also saves checkpoint data to `SimulationName.bak`. This is more than just current positions: the timestepping scheme needs (e.g.) the last two data points to interpolate forwards.
* To restart simulation from checkpoint, `./MPIFilament SimulationName` (where `SimulationName.bak` is the backup file).
* If `SimulationName` in `config.hpp` is the same as the checkpoint file, then the program will add data to the bottom of the existing `SimulationName.dat` data file.
* To start a simulation from time 0, using checkpoint backup data as position initialisation, use `./MPIFilament SimulationName restart`

### Restrictions ###
* Can not cope with filaments of different lengths.
* Can not cope with filaments with different swim phases.
