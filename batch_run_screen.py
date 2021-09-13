# Use this to run a load of simulations on this machine using screen.
#
# -- akt12  26/07/2019

import os
import sys
import fileinput
import subprocess
import datetime


def plural(word, num):
    if num == 1:
        return word
    else:
        return word + "s"

# 0 Self, 1 Timestamp, 2 Nodes, 3 Restart

if len(sys.argv) < 2:
    print("ERROR: You need to provide the batch timestamp as the first argument.")
    die()

if len(sys.argv) < 3:
    num_nodes = 4
else:
    num_nodes = int(sys.argv[2])

if len(sys.argv) == 4:
    force_restart_from_time_zero = True
else:
    force_restart_from_time_zero = False

batch_timestamp = sys.argv[1]
os.chdir(batch_timestamp)
subdirectories = [name for name in os.listdir(".") if os.path.isdir(name)]
os.chdir('..')

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
    restart_from_time_zero = ""
    if len(backup_files) > 0:
        continue_from_file = backup_files[-1][:-4]
        if len(parameter_files) == 0:
            restart_from_time_zero = " restart"
        if force_restart_from_time_zero:
            restart_from_time_zero = " restart"
    else:
        continue_from_file = ""


    os.chdir('..')
    # RUN THE MOTHERFUCKER
    # .. on one processor.
    cmd = "screen -dm bash -c 'export OPENBLAS_NUM_THREADS=1; mpiexec -np " + str(num_nodes) + " ./MPIFilament " + continue_from_file + restart_from_time_zero + ";'"
    os.system(cmd)
    # Change dir
    os.chdir('../..')



    if len(continue_from_file) > 0:
        if len(restart_from_time_zero) == 0:
            print(str(n+1) + ". Submitted: " + sim + ", on " + str(num_nodes) + " nodes, continuing from " + continue_from_file)
        else:
            print(str(n+1) + ". Submitted: " + sim + ", on " + str(num_nodes) + " nodes, initialised from " + continue_from_file)
    else:
        print(str(n+1) + ". Submitted: " + sim + ", on " + str(num_nodes) + " nodes")


print('')
num_jobs = len(subdirectories)

print("Complete: " + str(num_jobs) + " " + plural("job",num_jobs) + " running on " + plural("screen",num_jobs) + ", which will disappear when finished.")
print('')
