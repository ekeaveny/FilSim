# Use this to compile a load of simulations to run on the cluster. Enter the
# parameters desired at the top.
# The files are compiled into [timestamp]/SimName0/MPIFilament
# Output to be placed in      [timestamp]/SimNamo0/output/
#
# No need to run anything else before running this file. It does it all.
#
# -- sfs14/akt12  31/05/2019

import os
import sys
import fileinput
import subprocess
import datetime

# Defaults
StepsPerPeriod = 0
RandomMessNumNetworkFilaments = [1]
Nsw = [1]
TimeSteps = 0
Kap_NetworkFilaments = [1]
plot_steps = 0
StepsPerSettlingTime = 0
Bnumber = [1]
Barrier_FS_over_2a = [15*3.14159265358979]
Barrier_distance_cap = 2.0
WeightPerLength = [[0, 0, 0]]
SpermCurvatureDecay = False
SaveExtendedSwimmerData = False
randSeed = [2]
TOL = 1.e-6

# Disabled
NetworkSpringConstantFactor = [0]  # k in the spring force k(x-L)
NetworkSpringNaturalLength = 0  # L in the spring force k(x-L)
NetworkSpringConnectionRadius = [0]  # Look this far away from end of filament to create a link with another filament end
ProbabilityOfCreatingANetworkSpringLink = 0  # Form a spring link with this probability
EnableSpringLinks = False
ConnectedNodesNumSwimmingFilaments = 0
ConnectedNodesNumNodes = 1
ConnectedNodesConnectionRadius = 0
ConnectedNodesMinNodeSpacing = 0
ConnectedNodesSpacingAwayFromNode = 0

# Set config data here
#NPTS = [[384 * 1, 384 * 1, 384 * 8], [384 * 8, 384 * 1, 384 * 1]]  # [NPTS_X, NPTS_Y, NPTS_Z]
#NPTS = [[128,128,256]]
#NPTS = [[384,384,384]]
NPTS = [[288,288,288]]
#NPTS = [[4*1, 8 * 1, 8 * 1]]  # [NPTS_X, NPTS_Y, NPTS_Z]
Nworm = 15
ImplEulerSteps = 1
initialFrictionSteps = 0

SedimentationProblem = False
SwimProblem = False
TwoParticleProblem = True

TestInitialisation1 = 'false'
TestInitialisation2 = 'false'
TestInitialisation3 = 'false'
TestInitialisation4 = 'false'
TwoParticleInitialisation = 'false'


# Swimming
if SwimProblem:
    SimName0 = 'Swim'  # 'ConvCheckBDFb'
    GTTallBoxInitialisation = 'false'
    RandomMessToSwimThroughInitialisation = 'true'
    RandomMessToSwimThroughInitialisationAllowOverlap = 'false'

    Nworm_swimmer = 15

    # Disabled connective stuff
    ConnectedNodesInitialisation = False
    ConnectedNodesInitialisationAllowOverlap = False
    ConnectedNodesNumSwimmingFilaments = 0
    ConnectedNodesNumNodes = 1 #
    L_worm = Nworm_swimmer*2.2
    ConnectedNodesConnectionRadius = 1.2*L_worm
    ConnectedNodesMinNodeSpacing = 2.2*ConnectedNodesNumNodes**(-1./3.)*L_worm
    ConnectedNodesSpacingAwayFromNode = 0
    ConnectedNodesConnectionProbability = 0
    EnableSpringLinks = False
    NetworkSpringConstant = 0#[0,0.01,0.1]  # k in the spring force k*(x-L)
    NetworkSpringNaturalLength = 0  # L in the spring force k*(x-L)

    # Freeze the swimmer
    FrozenSwimmer = 'false'
    # WeightPerLength = [[0, 0, -1]] # Comment out for normal use; only for frozen swimmer.

    # Work out speed you'd get without collisions, but use gait you get with collisions?
    GaitReset = 'true'

    initialFrictionSteps = 0 #0 # 1000

    StepsPerPeriod = 300
    RandomMessNumNetworkFilaments = [400]  # [180, 300, 420]
    Nsw = [1 + i for i in RandomMessNumNetworkFilaments]
    TimeSteps = 4 * StepsPerPeriod #22  # Total timesteps
    Kap_NetworkFilaments = [1]  # [3.6, 36, 3600]
    plot_steps = int(StepsPerPeriod / 50)  # Save/checkpoint every plot_steps timesteps
    SpermCurvatureDecay = True
    SaveExtendedSwimmerData = True
    #randSeed = [2,3,4,5,6,7,8,9,10,11]
    randSeed = [2]
    Barrier_FS_over_2a = [15*3.14159265358979] # !!!!!!!!!!!!!!!!!!!!!!!!!!!
    TOL = 1.e-6

# Sedimenting
if SedimentationProblem:
    SimName0 = 'Sed'  # 'ConvCheckBDFb'
    GTTallBoxInitialisation = 'true'
    RandomMessToSwimThroughInitialisation = 'false'

    StepsPerSettlingTime = 300
    Nsw = [1]#[800]
    TimeSteps = StepsPerSettlingTime/30 # 600 * StepsPerSettlingTime  # Total timesteps
    Bnumber = [1e3]
    plot_steps = int(StepsPerSettlingTime / 100)  # Save/checkpoint every plot_steps timesteps
    WeightPerLength = [[0, 0, -1]]
    #WeightPerLength = [[-1, 0, 0]]

# Two-particle mobility tests in a sea of filaments
if TwoParticleProblem:
    SimName0 = 'TwoParticle'  # 'ConvCheckBDFb'
    GTTallBoxInitialisation = 'false'
    RandomMessToSwimThroughInitialisation = 'false'
    RandomMessToSwimThroughInitialisationTwoParticle = 'true'
    RandomMessToSwimThroughInitialisationAllowOverlap = 'false'
    TwoParticleInitialisation = 'true' # Additional initialisation

    # Irrelevant connective stuff
    Nworm_swimmer = 15
    ConnectedNodesInitialisation = False
    ConnectedNodesInitialisationAllowOverlap = False
    ConnectedNodesNumSwimmingFilaments = 0
    ConnectedNodesNumNodes = 1 #
    L_worm = Nworm_swimmer*2.2
    ConnectedNodesConnectionRadius = 1.2*L_worm
    ConnectedNodesMinNodeSpacing = 2.2*ConnectedNodesNumNodes**(-1./3.)*L_worm
    ConnectedNodesSpacingAwayFromNode = 0
    ConnectedNodesConnectionProbability = 0
    EnableSpringLinks = False
    NetworkSpringConstant = 0#[0,0.01,0.1]  # k in the spring force k*(x-L)
    NetworkSpringNaturalLength = 0  # L in the spring force k*(x-L)

    # Freeze the swimmer
    FrozenSwimmer = 'false'
    # WeightPerLength = [[0, 0, -1]] # Comment out for normal use; only for frozen swimmer.

    # Work out speed you'd get without collisions, but use gait you get with collisions?
    GaitReset = 'false'

    initialFrictionSteps = 0 #0 # 1000

    StepsPerPeriod = 300
    RandomMessNumNetworkFilaments = [400]#[801,1201,2001]  # [180, 300, 420]
    Nsw = [2 + i for i in RandomMessNumNetworkFilaments]
    TimeSteps = 4 * StepsPerPeriod #22  # Total timesteps
    Kap_NetworkFilaments = [1]#,1000]  # [3.6, 36, 3600]
    plot_steps = int(StepsPerPeriod / 50)  # Save/checkpoint every plot_steps timesteps
    SpermCurvatureDecay = False
    SaveExtendedSwimmerData = False
    SaveTwoParticleData = True
    #randSeed = [2,3,4,5,6,7,8,9,10,11]
    randSeed = [2]
    Barrier_FS_over_2a = [15*3.14159265358979] # !!!!!!!!!!!!!!!!!!!!!!!!!!!
    TOL = 1.e-6

time = datetime.datetime.now().strftime("%y%m%d%H%M")
os.system('mkdir ' + time)

# backup existing config.hpp so that we can restore it later.
if os.path.isfile('config.hpp'):
    os.system("cp config.hpp config-backup.hpp")
    config_backup = True
else:
    config_backup = False

print("Compiling batch of simulations and placing them in /" + time + "/ ...\n")

for i, Nsw_i in enumerate(Nsw):
    RandomMessNumNetworkFilaments_i = RandomMessNumNetworkFilaments[i]
    for Kap_NetworkFilaments_i in Kap_NetworkFilaments:
        for k, Bnumber_i in enumerate(Bnumber):
            for j, NPTS_i in enumerate(NPTS):
                for randSeed_i in randSeed:
                    Barrier_FS_over_2a_i = Barrier_FS_over_2a[k]
                    if SedimentationProblem:
                        SimName = SimName0 + "-N" + str(Nsw_i) + "-B" + str(Bnumber_i) + "-TallBox" + "-NPTSx" + str(NPTS_i[0]) + "y" + str(NPTS_i[1]) + "z" + str(NPTS_i[2])
                    if SwimProblem:
                        SimName = SimName0 + "-Nnetfil" + str(RandomMessNumNetworkFilaments_i) + "-Kap" + str(Kap_NetworkFilaments_i) + "-NPTSx" + str(NPTS_i[0]) + "y" + str(NPTS_i[1]) + "z" + str(NPTS_i[2]) + "-randSeed" + str(randSeed_i)
                        #SimName = SimName0 + "-Nnetfil" + str(RandomMessNumNetworkFilaments_i) + "-frame0"
                    if TwoParticleProblem:
                        SimName = SimName0 + "-Nnetfil" + str(RandomMessNumNetworkFilaments_i) + "-Kap" + str(Kap_NetworkFilaments_i) + "-randSeed" + str(randSeed_i)
                    SimName = SimName.replace(".", "p")
                    continue_from_file = ''
                    #continue_from_file = 'FreeNetwork-Nnetfil420-Kap3p6-NPTS1536'
                    # continue_from_file = 'FreeNetwork-Nnetfil' + str(RandomMessNumNetworkFilaments_i) + '-frame0'  # Leave blank if desired
                    if len(continue_from_file) > 0:
                        print("Compiling " + SimName + ", which continues from " + continue_from_file + " ...")
                    else:
                        print("Compiling " + SimName + " ...")

                    NPTS_X_i = NPTS_i[0]
                    NPTS_Y_i = NPTS_i[1]
                    NPTS_Z_i = NPTS_i[2]
                    xxRandomMessBoxSize_x_i = NPTS_X_i / 3.298629
                    xxRandomMessBoxSize_y_i = NPTS_Y_i / 3.298629
                    xxRandomMessBoxSize_z_i = NPTS_Z_i / 3.298629
                    ConnectedNodesBoxSize_x_i = NPTS_X_i / 3.298629
                    ConnectedNodesBoxSize_y_i = NPTS_Y_i / 3.298629
                    ConnectedNodesBoxSize_z_i = NPTS_Z_i / 3.298629
                    xxWeightPerLength_i = WeightPerLength[j]

                    # change config.hpp
                    os.system("cp config_template.hpp config.hpp")

                    with open('config.hpp', 'r') as file:
                        filedata = file.read()

                        filedata = filedata.replace('xxSedimentationProblem', ['false', 'true'][int(SedimentationProblem)])
                        filedata = filedata.replace('xxSwimProblem', ['false', 'true'][int(SwimProblem)])
                        filedata = filedata.replace('xxTwoParticleProblem', ['false', 'true'][int(TwoParticleProblem)])
                        filedata = filedata.replace('xxGTTallBoxInitialisation', GTTallBoxInitialisation)
                        filedata = filedata.replace('xxRandomMessToSwimThroughInitialisationAllowOverlap', RandomMessToSwimThroughInitialisationAllowOverlap)
                        filedata = filedata.replace('xxRandomMessToSwimThroughInitialisationTwoParticle', RandomMessToSwimThroughInitialisationTwoParticle)
                        filedata = filedata.replace('xxRandomMessToSwimThroughInitialisation', RandomMessToSwimThroughInitialisation)
                        filedata = filedata.replace('xxFrozenSwimmer', FrozenSwimmer)
                        filedata = filedata.replace('xxGaitReset', GaitReset)
                        filedata = filedata.replace('xxTwoParticleInitialisation', TwoParticleInitialisation)
                        filedata = filedata.replace('xxConnectedNodesInitialisationAllowOverlap', ['false', 'true'][int(ConnectedNodesInitialisationAllowOverlap)])
                        filedata = filedata.replace('xxConnectedNodesInitialisation', ['false', 'true'][int(ConnectedNodesInitialisation)])
                        filedata = filedata.replace('xxTestInitialisation1', TestInitialisation1)
                        filedata = filedata.replace('xxTestInitialisation2', TestInitialisation2)
                        filedata = filedata.replace('xxTestInitialisation3', TestInitialisation3)
                        filedata = filedata.replace('xxTestInitialisation4', TestInitialisation4)

                        filedata = filedata.replace('xxSimName', SimName)
                        filedata = filedata.replace('xxNworm', str(Nworm))
                        filedata = filedata.replace('xxNsw', str(Nsw[i]))
                        filedata = filedata.replace('xxNPTS_X', str(NPTS_X_i) + "L")
                        filedata = filedata.replace('xxNPTS_Y', str(NPTS_Y_i) + "L")
                        filedata = filedata.replace('xxNPTS_Z', str(NPTS_Z_i) + "L")
                        filedata = filedata.replace('xxTimeSteps', str(TimeSteps))
                        filedata = filedata.replace('xxplot_steps', str(plot_steps))
                        filedata = filedata.replace('xxImplEulerSteps', str(ImplEulerSteps))
                        filedata = filedata.replace('xxinitialFrictionSteps', str(initialFrictionSteps))

                        filedata = filedata.replace('xxRandomMessBoxSize_x', str(xxRandomMessBoxSize_x_i))
                        filedata = filedata.replace('xxRandomMessBoxSize_y', str(xxRandomMessBoxSize_y_i))
                        filedata = filedata.replace('xxRandomMessBoxSize_z', str(xxRandomMessBoxSize_z_i))
                        filedata = filedata.replace('xxRandomMessNumNetworkFilaments', str(RandomMessNumNetworkFilaments_i))
                        filedata = filedata.replace('xxRandomMessNumSwimmingFilaments', '0')
                        filedata = filedata.replace('xxKap_NetworkFilaments', str(Kap_NetworkFilaments_i))
                        filedata = filedata.replace('xxStepsPerPeriod', str(StepsPerPeriod))
                        filedata = filedata.replace('xxSpermCurvatureDecay', ['false', 'true'][int(SpermCurvatureDecay)])
                        filedata = filedata.replace('xxSaveExtendedSwimmerData', ['false', 'true'][int(SaveExtendedSwimmerData)])

                        filedata = filedata.replace('xxConnectedNodesNumSwimmingFilaments', str(ConnectedNodesNumSwimmingFilaments))
                        filedata = filedata.replace('xyNworm_swimmer', str(Nworm_swimmer))
                        filedata = filedata.replace('xxConnectedNodesNumNodes', str(ConnectedNodesNumNodes))
                        filedata = filedata.replace('xxConnectedNodesBoxSize_x', str(ConnectedNodesBoxSize_x_i))
                        filedata = filedata.replace('xxConnectedNodesBoxSize_y', str(ConnectedNodesBoxSize_y_i))
                        filedata = filedata.replace('xxConnectedNodesBoxSize_z', str(ConnectedNodesBoxSize_z_i))
                        filedata = filedata.replace('xxConnectedNodesConnectionRadius', str(ConnectedNodesConnectionRadius))
                        filedata = filedata.replace('xxConnectedNodesMinNodeSpacing', str(ConnectedNodesMinNodeSpacing))
                        filedata = filedata.replace('xxConnectedNodesSpacingAwayFromNode', str(ConnectedNodesSpacingAwayFromNode))
                        filedata = filedata.replace('xxConnectedNodesConnectionProbability', str(ConnectedNodesConnectionProbability))

                        filedata = filedata.replace('xxWeightPerLength_x', str(xxWeightPerLength_i[0]))
                        filedata = filedata.replace('xxWeightPerLength_y', str(xxWeightPerLength_i[1]))
                        filedata = filedata.replace('xxWeightPerLength_z', str(xxWeightPerLength_i[2]))

                        filedata = filedata.replace('xxBnumber', str(Bnumber_i))
                        filedata = filedata.replace('xxBarrier_FS_over_2a', str(Barrier_FS_over_2a_i))
                        filedata = filedata.replace('xxBarrier_distance_cap', str(Barrier_distance_cap))
                        filedata = filedata.replace('xxStepsPerSettlingTime', str(StepsPerSettlingTime))

                        # Network springs
                        filedata = filedata.replace('xxEnableSpringLinks', ['false', 'true'][int(EnableSpringLinks)])
                        filedata = filedata.replace('xxNetworkSpringConstant', str(NetworkSpringConstant))
                        filedata = filedata.replace('xxNetworkSpringNaturalLength', str(NetworkSpringNaturalLength))

                        filedata = filedata.replace('xxrandSeed', str(randSeed_i))

                        filedata = filedata.replace('xxTOL', str(TOL))

                    with open('config.hpp', 'w') as file:
                        file.write(filedata)

                    # make
                    #os.system("export OPENBLAS_NUM_THREADS=1")  # sets number of threads openBLAS uses to 1 (no communication overhead between cpus)
                    #subprocess.call("./make_wrapper_cx2.sh")  # Make
                    #process = subprocess.Popen(["make","MPI_FIL_debug"], stdout=subprocess.PIPE)
                    process = subprocess.Popen(["make","MPI_FIL_twoparticle"], stdout=subprocess.PIPE)
                    #print "NOTE DEBUG MODE IS ON"
                    process.wait()

                    # make new subdirectory called SimName, copy binary file
                    os.system('mkdir ' + time + "/" + SimName + ' && cp MPIFilament ' + time + "/" + SimName)
                    os.chdir(time + "/" + SimName)
                    os.system('mkdir output')
                    os.chdir("../..")
                    os.system("cp config.hpp " + time + "/" + SimName)
                    if len(continue_from_file) > 0:
                        os.system("cp ./output/" + continue_from_file + ".par " + time + "/" + SimName + "/output")
                        os.system("cp ./output/" + continue_from_file + ".dat " + time + "/" + SimName + "/output")
                        os.system("cp ./output/" + continue_from_file + ".bak " + time + "/" + SimName + "/output")

# Return backup config.hpp
if config_backup:
    os.system("cp config-backup.hpp config.hpp")
    os.system("rm config-backup.hpp")

print("\nComplete. Simulations have been compiled and placed in /" + time + "/")
#print("\nUse command `python batch_run.py " + time + "" + '` to run.')
print("\nUse command `python batch_run_screen.py " + time + "" + '` to run.')
