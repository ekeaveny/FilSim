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

COLOR_GREEN = u"\u001b[42m\u001b[30m"
COLOR_YELLOW = u"\u001b[43m\u001b[30m"
COLOR_RED = u"\u001b[41m"
COLOR_RESET = u"\u001b[0m\u001b[0m"

batch_timestamp = sys.argv[1]
os.chdir(batch_timestamp)
subdirectories = [name for name in os.listdir(".") if os.path.isdir(name)]
os.chdir('..')

subdirectories.sort()
num_finished_sims = 0
num_sims_in_progress = 0
# print('------------------------------------------------------------------')
for n, sim in enumerate(subdirectories):
    # Find name of continuation file, if it exists
    os.chdir(batch_timestamp + "/" + sim + "/output")
    files = [name for name in os.listdir(".") if os.path.isfile(name)]
    # Sort by modified time, so we pick the last edited .bak file
    files.sort(key=lambda x: os.path.getmtime(x))
    backup_files = [name for name in files if name[-4:] == '.bak']
    parameter_files = [name for name in files if name[-4:] == '.par']
    data_files = [name for name in files if name[-4:]
                  == '.dat' and name[-16:] != "-springlinks.dat"]
    data_files_modified_time = [os.path.getmtime(x) for x in data_files]
    restart_from_time_zero = ""
    finished = False
    in_progress = False

    if len(data_files) > 0:
        # Check end of data file to see whether it is complete
        no_data_files = False
        try:
            PY3 = sys.version_info[0] == 3
            data_file = data_files[0]
            if PY3:
                last_nt = int(subprocess.check_output(
                    ['tail', '-1', data_file], encoding='utf8').split(" ")[0])  # CX2 version
                #last_nt = int(subprocess.check_output(['tail', '-1', '-r', data_file], encoding='utf8').split(" ")[0])
            else:
                last_nt = int(subprocess.check_output(
                    ['tail', '-1', data_file]).split(" ")[0])  # CX2 version
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
                num_finished_sims += 1
            else:
                total_len = len(str(requested_timesteps))
                progress = str(last_nt).rjust(total_len, ' ') + \
                    "/" + str(requested_timesteps)
            modified_time = data_files_modified_time[0]
            thirty_minutes = 30 * 60
            if modified_time >= time.time() - thirty_minutes:
                in_progress = True
                num_sims_in_progress += 1
        except:
            print("There was a problem testing whether simulation " + str(n) + " has already finished.")
            progress = "-ERROR-"
    else:
        no_data_files = True

    os.chdir('../../..')

    if finished:
        flag = COLOR_GREEN + " FINISHED  " + COLOR_RESET
    elif no_data_files:
        flag = COLOR_RED + "NOT STARTED" + COLOR_RESET
    else:
        flag = COLOR_YELLOW + progress + COLOR_RESET

    if in_progress:
        prog_flag = " " + COLOR_YELLOW + "!" + COLOR_RESET + " "
    else:
        prog_flag = "   "

    print flag + prog_flag + sim

if num_sims_in_progress > 0:
    print(COLOR_YELLOW + "!" + COLOR_RESET + ' = Data file written to in last 30 mins')
else:
    print("No data files written to in the last 30 mins.")
