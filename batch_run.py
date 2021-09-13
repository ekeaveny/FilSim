# Use this to run a load of simulations to run on the cluster.
# Expects that you've compiled the simulations using  batch_compile.py .
#
# Will continue from whatever .bak file is in the output directory
#
# No need to do anything else before running this file. It does it all.
#
# WARNING: This will overwrite  pbs_job_script_run.sh
#
# -- akt12  28/05/2019

import os
import sys
import fileinput
import subprocess
import datetime
import time


def plural(word, num):
    if num == 1:
        return word
    else:
        return word + "s"


if len(sys.argv) < 2:
    print("ERROR: You need to provide the batch timestamp as the first argument.")
    sys.exit()

if len(sys.argv) == 3:
    force_restart_from_time_zero = True
else:
    force_restart_from_time_zero = False

debug_mode_do_not_submit_anything = False

do_not_submit_if_data_file_modified_less_than_num_secs_ago = 30*60

# batch_walltime = '47:59:00'
# batch_walltime = '0:59:00'
#batch_walltime = '23:59:00' # normal
#batch_walltime = '29:59:00' # normal
# batch_walltime = '1:59:00'
batch_walltime = '0:05:00'
#batch_walltime = '1:00:00'
#batch_num_nodes = '64'  # Number of nodes (individual machines)
#batch_num_nodes = '2'  # Number of nodes (individual machines) (normal)
batch_num_nodes = '6'  # Number of nodes (individual machines)
#batch_walltime = '1:59:00'
#batch_num_nodes = '3'
batch_num_cpus_per_node = '24'  # Num CPUs/node to use. Must be 24 for cx2.
# Memory upper bound: NPTSxy*NPTSxy*NPTSz*10*8/10^9/batch_num_nodes (in GB)
#batch_memory = '120gb'  # Memory requested per node. Up to 120gb.
# batch_memory_num = 150 # 16GB/node
batch_memory_num = 32 # 16GB/node
#batch_memory_num = 16 # 16GB/node <----- normal choice
#batch_mpiprocs = '144'
express = False
'''
express = True
# Num CPUs/node to use. Must be 256 for express.
batch_num_cpus_per_node = '256'
batch_num_nodes = '1'  # Number of nodes (individual machines)
batch_memory_num = 960 # 16GB/node
'''
batch_mpiprocs = batch_num_cpus_per_node  # By default.

if debug_mode_do_not_submit_anything:
    print('!!! DEBUG MODE - NOTHING WILL BE SUBMITTED !!!')

batch_timestamp = sys.argv[1]
os.chdir(batch_timestamp)
subdirectories = [name for name in os.listdir(".") if os.path.isdir(name)]
os.chdir('..')

num_finished_sims = 0
num_sims_in_progress = 0
print('')
#print('------------------------------------------------------------------')
for n, sim in enumerate(subdirectories):
    # Find name of continuation file, if it exists
    os.chdir(batch_timestamp + "/" + sim + "/output")
    files = [name for name in os.listdir(".") if os.path.isfile(name)]
    # Sort by modified time, so we pick the last edited .bak file
    files.sort(key=lambda x: os.path.getmtime(x))
    backup_files = [name for name in files if name[-4:] == '.bak']
    parameter_files = [name for name in files if name[-4:] == '.par']
    data_files = [name for name in files if name[-4:] == '.dat' and name[-16:] != "-springlinks.dat"]
    data_files_modified_time = [os.path.getmtime(x) for x in data_files]
    restart_from_time_zero = ""
    finished = False
    in_progress = False
    if len(backup_files) > 0:
        continue_from_file = backup_files[-1][:-4]
        if len(parameter_files) == 0:
            restart_from_time_zero = "restart"
        if len(data_files) == 0:
            restart_from_time_zero = "restart"
        if force_restart_from_time_zero:
            restart_from_time_zero = " restart"
        elif len(data_files) > 0:
            # Check end of data file to see whether it is complete
            try:
                PY3 = sys.version_info[0] == 3
                data_file = data_files[0]
                if PY3:
                    last_nt = int(subprocess.check_output(['tail', '-1', data_file], encoding='utf8').split(" ")[0]) # CX2 version
                    #last_nt = int(subprocess.check_output(['tail', '-1', '-r', data_file], encoding='utf8').split(" ")[0])
                else:
                    last_nt = int(subprocess.check_output(['tail', '-1', data_file]).split(" ")[0]) # CX2 version
                    #last_nt = int(subprocess.check_output(['tail', '-1', '-r', data_file]).split(" ")[0])
                # Find number of timesteps asked for
                with open(parameter_files[0]) as fp:
                    for i, line in enumerate(fp):
                        if i == 1:
                            data_parameters = line.split(" ")
                        elif i > 1:
                            break
                requested_timesteps = int(data_parameters[13])
                if last_nt == requested_timesteps:
                    finished = True
                    num_finished_sims+=1

                if not finished:
                    modified_time = data_files_modified_time[0]
                    if modified_time >= time.time() - do_not_submit_if_data_file_modified_less_than_num_secs_ago:
                        in_progress = True
                        num_sims_in_progress += 1
            except:
                print("There was a problem testing whether simulation " + str(n) + " has already finished or is currently running.")


    else:
        continue_from_file = ""

    os.chdir('../../..')

    # Find config.hpp number of grid points so we know how much memory to ask
    # for.
    os.chdir(batch_timestamp + "/" + sim)
    npts_x = 0
    npts_y = 0
    npts_z = 0
    with open("config.hpp", "r") as configfile:
        for line in configfile:
            if line[0:14] == "#define NPTS_X":
                npts_x = int(line[15:25].split("L")[0])
            if line[0:14] == "#define NPTS_Y":
                npts_y = int(line[15:25].split("L")[0])
            if line[0:14] == "#define NPTS_Z":
                npts_z = int(line[15:25].split("L")[0])
    os.chdir('../..')
    if npts_x > 0 and npts_y > 0 and npts_z > 0:
        #batch_memory_num = min(8,int(npts_x*npts_y*npts_z*10*8*6*2*5/1e9/int(batch_num_nodes))+1)
        #batch_memory_num = 16
        batch_memory = str(batch_memory_num) + "gb"
        # the 6 is magic
        # the 2*5 is for good luck

        if (npts_x*npts_y*npts_z) % (int(batch_num_cpus_per_node)*int(batch_num_nodes)) != 0:
            error_message = "NPTS_X*NPTS_Y*NPTS_Z must be divisible by the "
            error_message += "number of nodes you are requesting. "
            error_message += "NPTS_X = " + str(npts_x) + ", "
            error_message += "NPTS_Y = " + str(npts_y) + ", "
            error_message += "NPTS_Z = " + str(npts_z) + " and you have requested "
            error_message += str(int(batch_num_cpus_per_node)*int(batch_num_nodes)) + " nodes."
            sys.exit(error_message)

    os.system("cp pbs_job_script_run_template.sh pbs_job_script_run.sh")
    with open('pbs_job_script_run.sh', 'r') as file:
        filedata = file.read()

        filedata = filedata.replace('xxWALLTIME', batch_walltime)
        filedata = filedata.replace('xxNUMNODES', batch_num_nodes)
        filedata = filedata.replace('xxNUMCPUS', batch_num_cpus_per_node)
        filedata = filedata.replace('xxMPIPROCS', batch_mpiprocs)
        filedata = filedata.replace('xxMEMORY', batch_memory)
        filedata = filedata.replace('xxPWD', os.getcwd())
        filedata = filedata.replace('xxTIMESTAMP', batch_timestamp)
        filedata = filedata.replace('xxSIM', sim)
        filedata = filedata.replace('xxCONTINUEFROM', continue_from_file)
        filedata = filedata.replace('xxRESTARTFROMTIMEZERO', restart_from_time_zero)

    with open('pbs_job_script_run.sh', 'w') as file:
        file.write(filedata)

    if finished:
        print(str(n+1) + ". Skipped:   " + sim + " is already complete.")
    elif in_progress:
        print(str(n+1) + ". Skipped:   " + sim + " is currently running.")
    else:
        if len(continue_from_file) > 0:
            if len(restart_from_time_zero) == 0:
                print(str(n+1) + ". Submitted: " + sim + ", continuing from " + continue_from_file)
            else:
                print(str(n+1) + ". Submitted: " + sim + ", initialised from " + continue_from_file)
        else:
            print(str(n+1) + ". Submitted: " + sim)

        print("   Memory requested: " + batch_num_nodes + "x" + batch_memory)

        if not debug_mode_do_not_submit_anything:
            if express:
                os.system('qsub -q express -P exp-00045 pbs_job_script_run.sh')
            else:
                os.system('qsub pbs_job_script_run.sh ')


print('')
num_jobs = len(subdirectories)

num_submitted_jobs = num_jobs-num_finished_sims-num_sims_in_progress
print("Complete: " + str(num_submitted_jobs) + " " + plural("job",num_submitted_jobs) + " submitted to " + batch_num_nodes + "x" + batch_num_cpus_per_node + " " + plural("processor",int(batch_num_nodes)*int(batch_num_cpus_per_node)) + " with walltime of " + batch_walltime + ".")
if num_finished_sims > 0:
    print("          " + str(num_finished_sims) + " " + plural("job",num_finished_sims) + " skipped because they have already finished.")
if num_sims_in_progress > 0:
    print("          " + str(num_sims_in_progress) + " " + plural("job",num_sims_in_progress) + " skipped because their data files were written to in the last " + str(int(do_not_submit_if_data_file_modified_less_than_num_secs_ago/60)) + " mins.")
#print('------------------------------------------------------------------')
print('')
print('To see queue:                      `qstat`')
#print('To see queue with estd start time: `qstat -w -T`')
#print('To see queue with details in text: `qstat -f`')
print('To see availability of resources:  `availability`')
print('To delete all jobs by me:          `qselect -u akt12 | xargs qdel`')
#print('------------------------------------------------------------------')
print('')
if debug_mode_do_not_submit_anything:
    print('!!! DEBUG MODE - NOTHING HAS ACTUALLY BEEN SUBMITTED !!!')
