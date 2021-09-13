"""
Plots swimmer speeds from C++ simulations.

Looks through output folders, sees which have the same configs, and then
averages them if appropriate.

AKT 15/04/2019. (Python 3)
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess
import os

output_folders = ['../output/1912102311/',
                  '../output/1912191427/',
                  '../output/1912232152/',
                  '../output/2001061917/',
                  '../output/2001201558/',
                  '../output/2001281337/',
                  '../output/2002141041/',
                  '../output/2002170830/',
                  '../output/2003091939/'
                  ]

# Find all simulations in the above folders
all_sim_home_paths = []
all_sim_names = []
for folder in output_folders:
    os.chdir(folder)
    subdirectories = [name for name in os.listdir(".") if os.path.isdir(name)] # SimName
    for sim in subdirectories:
        all_sim_home_paths.append(folder+sim) # ../output/1910blah/SimName
        all_sim_names.append(sim) # SimName
    os.chdir('../../postprocessing')

# Find parameters and see which should be averaged over
parameters = []
for n, name in enumerate(all_sim_names):
    output_folder = all_sim_home_paths[n] + '/output/'
    parameter_file = output_folder + name + '.par'

    # Read in parameters
    with open(parameter_file) as fp:
        for i, line in enumerate(fp):
            if i == 1:
                data_parameters = line.split(" ")
            elif i > 1:
                break

    N_sw = int(data_parameters[2])
    N_w = int(data_parameters[3])
    #Kap_NetworkFilaments = float(data_parameters[24])
    Kap_NetworkFilaments = float(name.split("-")[2][3:])
    if len(name.split("-")) > 4:
        randSeed = int(name.split("-")[4][8:])
    else:
        randSeed = 2
    #LfcmBox_x = float(data_parameters[19])
    #LfcmBox_y = float(data_parameters[25])
    #LfcmBox_z = float(data_parameters[24])
    LfcmBox_x = float(data_parameters[19])
    LfcmBox_y = LfcmBox_x
    LfcmBox_z = LfcmBox_x # <-------------------------------------------------------APPEARS TWICE
    phi = 4./3.*3.1415926*(N_sw-1)*N_w/(LfcmBox_x*LfcmBox_y*LfcmBox_z)*100 # as percentage

    parameters.append([N_sw,Kap_NetworkFilaments,randSeed,phi,n])

    if n == 0:
        # Titles
        print ("     [" + ("Name" + " "*30)[0:30] + "][" +
                    ("M" + " "*4)[0:4] + "][" +
                    ("Kbnet" + " "*5)[0:5] + "][" +
                    ("randSeed" + " "*8)[0:8] + "]")
    print (     "%4i" % n + " " +
                "[" + (output_folder + name)[0:30] + "][" +
                "%4i" % N_sw + "][" +
                "%5i" % Kap_NetworkFilaments + "][" +
                "%8i" % randSeed + "]")


parameters.sort()
order = [i[-1] for i in parameters]

#from IPython import embed
#embed()

swimmer_velocity = []
timess = []

for n in order:
    name = all_sim_names[n]
    output_folder = all_sim_home_paths[n] + '/output/'

    parameter_file = output_folder + name + '.par'
    #data_file = output_folder + name + '-every10throw.dat' # '.dat'#
    data_file = output_folder + name + '.dat'#
    swimmer_file = output_folder + name + '-swimmer.dat'#

    plot_every_n_frames = 10

    # Read in parameters
    with open(parameter_file) as fp:
        for i, line in enumerate(fp):
            if i == 1:
                data_parameters = line.split(" ")
            elif i > 1:
                break

    # Read in data
    if name[-4:] == ".dat":
        name = name[:-4]
    print("Reading in from " + output_folder + name + '.dat' + ' ...')

    N_sw = int(data_parameters[2])
    N_w = int(data_parameters[3])
    Sp4 = float(data_parameters[7])
    KAP = float(data_parameters[8])
    #Kap_NetworkFilaments = float(data_parameters[24])
    Kap_NetworkFilaments = float(name.split("-")[2][3:])
    mu = float(data_parameters[10])
    L = int(data_parameters[3])*float(data_parameters[6])
    StepsPerPeriod = int(data_parameters[14]) # 300
    print StepsPerPeriod
    period = StepsPerPeriod
    omega = (KAP*(Sp4)/(4*3.1415926*L**4));
    f = omega/(2*3.1415926);
    dt = 1/(period*f);
    T = dt*period
    B = 1000
    LfcmBox_x = float(data_parameters[19])
    LfcmBox_y = LfcmBox_x
    LfcmBox_z = LfcmBox_x # <-------------------------------------------------------APPEARS TWICE
    #LfcmBox_y = float(data_parameters[25])
    #LfcmBox_z = float(data_parameters[24])
    phi = 4./3.*3.1415926*(N_sw-1)*N_w/(LfcmBox_x*LfcmBox_y*LfcmBox_z)*100 # as percentage
    #B = int(data_parameters[26])
    if int(LfcmBox_z) == 0:
        LfcmBox_z = LfcmBox_x
    if int(LfcmBox_y) == 0:
        LfcmBox_y = LfcmBox_x
    #data_saved_every_n_timesteps = int(data_parameters[20]) *10# <===========================================
    data_saved_every_n_timesteps = 1 # HARDCODED to 1 for swimmer data # int(data_parameters[20]) # <===========================================

    final_nt = int(subprocess.check_output(['tail', '-1', '-r', swimmer_file]).split(" ")[0])
    first_nt = 0  # Could replace with head if necessary
    num_rows = (final_nt - first_nt) // data_saved_every_n_timesteps + 1

    num_rows = num_rows -1#<---------------------------------------------

    num_frames = num_rows // plot_every_n_frames

    print "Number of rows of data:       ", num_rows
    print "Data saved every n timesteps: ", data_saved_every_n_timesteps
    print "Number of frames to plot:     ", num_frames

    times = []
    Us = []
    UUs = []
    com_filament_1 = []
    nts = []
    U_coms = []
    U_dot_ts = []
    t_avgs = []
    UUUs = []
    UUUUs = []
    ts = []

    with open(swimmer_file, "r") as datafile:
        file_row_number = -1  # Because we want to ignore header row, right?
        previous_row_num = -1
        previous_nt = -1

        for frame_num in range(0,num_frames):
            row_num = int(frame_num * plot_every_n_frames)
            if row_num != previous_row_num:  # Frame 0 is generated twice by the video generator. Annoying, isn't it?
                while file_row_number < row_num:
                    datafile.readline()
                    file_row_number += 1
                row = datafile.readline()
                file_row_number += 1
                row = row.split(" ")
            previous_row_num = row_num

            try:
            #if 1==1:
                nt = int(row[0])

                print "frame is", frame_num, ". row_num is", row_num, ". nt is", nt
                print LfcmBox_x/L

                X = np.array(row[1:7 * N_w + 1:7], dtype=np.float32)
                Y = np.array(row[2:7 * N_w + 2:7], dtype=np.float32)
                Z = np.array(row[3:7 * N_w + 3:7], dtype=np.float32)

                q0 = np.array(row[4:7 * N_w + 4:7], dtype=np.float32)
                q1 = np.array(row[5:7 * N_w + 5:7], dtype=np.float32)
                q2 = np.array(row[6:7 * N_w + 6:7], dtype=np.float32)
                q3 = np.array(row[7:7 * N_w + 7:7], dtype=np.float32)

                Ux = np.array(row[7 * N_w + 1:16 * N_w + 1:9], dtype=np.float32)
                Uy = np.array(row[7 * N_w + 2:16 * N_w + 2:9], dtype=np.float32)
                Uz = np.array(row[7 * N_w + 3:16 * N_w + 3:9], dtype=np.float32)

                X_com = np.array([X.mean(), Y.mean(), Z.mean()])

                U_com = np.array([Ux.mean(), Uy.mean(), Uz.mean()])
                t1 = 1 - 2*q2**2 - 2*q3**2
                t2 = 2*(q1*q2 + q3*q0)
                t3 = 2*(q1*q3 - q2*q0)
                t_avg = np.array([t1.mean(), t2.mean(), t3.mean()])
                t_avg_hat = -t_avg/np.linalg.norm(t_avg) # Note t points towards tail

                t = nt*dt/T

                if nt > 0:
                    U_dot_ts.append(np.dot(U_com,t_avg_hat))
                    U_coms.append(np.linalg.norm(U_com))
                    com_filament_1.append(X_com)
                    nts.append(nt)
                    ts.append(t)
                    UUUs.append(U_com)
                    t_avgs.append(t_avg_hat)



                period_in_frames = period//plot_every_n_frames//data_saved_every_n_timesteps
                one_period_ago = -1-period_in_frames
                xxx = 0


                if len(nts) >= -one_period_ago:

                    avg_t_over_last_period = np.array(t_avgs[one_period_ago:-1]).mean(axis=0)

                    UUU = np.dot(np.array(UUUs[one_period_ago:-1]),avg_t_over_last_period).mean()
                    UUUUs.append(UUU/(L/T))

                    print "time: ", t, ", speed: ", UUU, ", dt/T: ", dt/T

                    times.append(t)

                previous_nt = nt
            except:
                print "ERRO"



    timess.append(times)
    swimmer_velocity.append(UUUUs)

#from IPython import embed
#embed()

previous_parameters = []
parameters_unique = []
for p, para in enumerate(parameters):
    if para[:2] == previous_parameters:
        parameters_unique[-1].append(p)
    else:
        parameters_unique.append([p])
    previous_parameters = para[:2]

for n, p in enumerate(parameters_unique[1:]):
    N_sw = parameters[p[0]][0]
    Kap_NetworkFilaments = parameters[p[0]][1]

    #if N_sw >= 2000:
    #    from IPython import embed
    #    embed()

    min_length = min([len(timess[q]) for q in p])
    ti = timess[p[0]][:min_length]
    unhindered_swimmer_constant_speed = np.mean(swimmer_velocity[0][300:])
    #from IPython import embed
    #embed()
    u = np.mean(np.array([np.array(swimmer_velocity[q][:min_length]) for q in p]),axis=0) / unhindered_swimmer_constant_speed

    label = "$\\phi=" + "%.1f"%parameters[p[0]][3] + "\%$, $K_B=" + "%i"%parameters[p[0]][1] + "$ ($n=" + str(len(p)) + "$)"
    color='C' + str(n//2)
    dashes = [(1,0),(1,1)][n%2]
    plt.plot(ti,u,label=label,color=color,dashes=dashes)


#        times[pp]
#        swimmer_velocity[pp]
    #swimmer_velocity[0]

'''
for n in range(len(parameters)-1):
    if parameters[n][:-1] ==

if n == 0:
    UUUUs_ref = UUUUs
else:
    n = n - 1
    u = [UUUUs[i]/UUUUs_ref[i] for i in range(len(UUUUs))]
    color='C' + str(n//2)
    dashes = [(1,0),(1,1)][n%2]
    plt.plot(times,u,label=label,color=color,dashes=dashes)
'''
time_taken_for_unhindered_swimmer_to_swim_box_length = LfcmBox_x/L/unhindered_swimmer_constant_speed
for i in range(3):
    plt.axvline(x=i*time_taken_for_unhindered_swimmer_to_swim_box_length)
plt.xlim([0,50])
plt.xlabel('time/T')
plt.ylabel('swimmer speed/speed in empty fluid')
plt.grid()
plt.legend(loc='best')
plt.show()
