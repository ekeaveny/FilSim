// Template for config.hpp to use with batch_compile.py

// Set all simulation paramters

// inclusion guards
#ifndef CONFIG_IS_INCLUDED
#define CONFIG_IS_INCLUDED

// for DEBUGGING [if all false, full functionality guaranteed]
#define verbose false
#define MPIverbose false
#define flush_on false // Dynamic newlines: good for screen viewing, not saving.
#define notAllTorques false
#define PROFILING true

// #define myOpenMPthreads 4 // For RPY only.

// Sedimentation problem or swimming problem for unit choice
#define SwimProblem xxSwimProblem
#define SedimentationProblem xxSedimentationProblem
#define TwoParticleProblem xxTwoParticleProblem

// Initial setup
#define GTinitialisation false
#define RANDinitialisation false
#define StackedInitialisation false
#define GTTallBoxInitialisation xxGTTallBoxInitialisation

#define RandomMessToSwimThroughInitialisation xxRandomMessToSwimThroughInitialisation
#define RandomMessToSwimThroughInitialisationTwoParticle xxRandomMessToSwimThroughInitialisationTwoParticle
#define RandomMessToSwimThroughInitialisationAllowOverlap xxRandomMessToSwimThroughInitialisationAllowOverlap
#define FrozenSwimmer xxFrozenSwimmer
#define GaitReset xxGaitReset

#define TwoParticleInitialisation xxTwoParticleInitialisation
#define TwoParticleSeparation 0

#define RandomMessNumNetworkFilaments xxRandomMessNumNetworkFilaments
#define RandomMessNumSwimmingFilaments xxRandomMessNumSwimmingFilaments
// Will put centre of filament inside a box of size
// [0,RandomMessBoxSize_xy]*[0,RandomMessBoxSize_xy]*[0,RandomMessBoxSize_z]
#define RandomMessBoxSize_x xxRandomMessBoxSize_x // 465.648 // 38.8 // 77.60
#define RandomMessBoxSize_y xxRandomMessBoxSize_y // 465.648 // 38.8 // 77.60
                               // 77.608 equivalent to box with NPTS=256. 116.412 equivalent to NPTS=384
#define RandomMessBoxSize_z xxRandomMessBoxSize_z

#define ConnectedNodesInitialisation xxConnectedNodesInitialisation
#define ConnectedNodesInitialisationAllowOverlap xxConnectedNodesInitialisationAllowOverlap

#define TestInitialisation1 xxTestInitialisation1
#define TestInitialisation2 xxTestInitialisation2
#define TestInitialisation3 xxTestInitialisation3
#define TestInitialisation4 xxTestInitialisation4


// Number of filaments/swimmers
// in the present implementation, this number is limited by the stack size, if ~STACKALLOCATION
// For lattices, note that for Nx,Ny,Nz grid points, num of filaments required = 3 Nx Ny Nz - Nx Ny - Ny Nz - Nx Nz.
//#define Nsw (3*LatticeNi*LatticeNj*LatticeNk - LatticeNi*LatticeNj - LatticeNj*LatticeNk - LatticeNj*LatticeNi + LatticeNumFilamentsOutside) // = 54 (for lattice) + 1
#define Nsw xxNsw // (RandomMessNumNetworkFilaments+RandomMessNumSwimmingFilaments)

// Number of particles per filament/swimmer
#define Nworm xxNworm // Default filament length
#define max_Nworm xxNworm // Maximum possible filament length (for array allocation)
#define Nworm_swimmer xyNworm_swimmer // Assume all swimming filaments are the same length. If you want to change this, be careful of the definition of omega

// Number of filaments/swimmers
//#define Nsw 10
// Number of particles per filament/swimmer
//#define Nworm 20

// FCM configuration
#define NPTS_X xxNPTS_X // 128 // 256
#define NPTS_Y xxNPTS_Y // 128 // 256
#define NPTS_Z xxNPTS_Z // 128 // 256
#define ngd 20//15

// constants
#define a  1.0 // segment radius
#define ah 1.0 // head radius
#define DL 2.2 // segment distance
#define MU 1.0 // viscosity
#define WeightPerLengthX xxWeightPerLength_x // gravitational force line density in +X-direction
#define WeightPerLengthY xxWeightPerLength_y // gravitational force line density in +Y-direction
#define WeightPerLengthZ xxWeightPerLength_z // gravitational force line density in +Z-direction
// #define Lbox 2. // initial box size in swimmer lengths
#define InitialSeparation 50. //separation in GT initialisation

// Random seed
#define randSeed xxrandSeed

// Sperm number and bending modulus
#define Sp4       1187 // now using def of Sp = L(4piomegamu/Kb)^0.25
//           24.0 // using other def of Sp involving 2pi // not Sp, but Sp^4
#define Bnumber   xxBnumber  //1.e-3 // elasto-gravitational number
#define Barrier_FS_over_2a   xxBarrier_FS_over_2a
#define Barrier_distance_cap   xxBarrier_distance_cap

#if SwimProblem || TwoParticleProblem
    #define Kap                   1000    //1800.0 // in terms of Bnumber, Kap = WL^3/B
    #define C                     Kap     //1800.0
    #define Kap_NetworkFilaments  xxKap_NetworkFilaments      // (Kap=36 corresponds to B=1000. Kap=36000 corresponds to B=1. Of course B is not really well defined without gravity.)
    #define C_NetworkFilaments    Kap_NetworkFilaments
#endif

// Target curvature per length (double!)
#define curvature 8.25

// Swimming: Wave number, helical swimming alpha and beta
#define K0 0.75 // 1
#define SwimmingHelixAlpha 0.
#define SwimmingHelixBeta 1.
#define SpermCurvatureDecay xxSpermCurvatureDecay // Swims like a sperm, with decaying preferred curvature after half length

// Networks
#define EnableSpringLinks xxEnableSpringLinks
#define NetworkSpringConstant xxNetworkSpringConstant // k in the spring force k(x-L)
#define NetworkSpringNaturalLength xxNetworkSpringNaturalLength // L in the spring force k(x-L)

// Connected nodes
#define ConnectedNodesNumNodes xxConnectedNodesNumNodes
#define ConnectedNodesConnectionRadius xxConnectedNodesConnectionRadius
#define ConnectedNodesMinNodeSpacing xxConnectedNodesMinNodeSpacing
#define ConnectedNodesSpacingAwayFromNode xxConnectedNodesSpacingAwayFromNode
#define ConnectedNodesConnectionProbability xxConnectedNodesConnectionProbability

#define ConnectedNodesNumSwimmingFilaments xxConnectedNodesNumSwimmingFilaments
#define ConnectedNodesBoxSize_x xxConnectedNodesBoxSize_x
#define ConnectedNodesBoxSize_y xxConnectedNodesBoxSize_y
#define ConnectedNodesBoxSize_z xxConnectedNodesBoxSize_z

// Time integration
#define StepsPerPeriod xxStepsPerPeriod  // for swimmer problem
                            // note StepsPerPeriod should be a multiple of 100 unless you redefine plot_steps below.
#define StepsPerSettlingTime xxStepsPerSettlingTime // for sedimentation problem
#define TimeSteps xxTimeSteps//      StepsPerPeriod/30 // 30000

// Strain twist vectors
#define STRAINTWISTX 0.
#define STRAINTWISTY 0.
#define STRAINTWISTZ 0.

// number of initial time steps computed with implicit euler instead of BDF scheme
#define ImplEulerSteps xxImplEulerSteps//2//100000 //1
#define initialFrictionSteps xxinitialFrictionSteps // number of steps with local friction only

// Tolerance for non-linear system, gmres, and maximum number of iterations
#define TOL xxTOL //10 //1.e-8 //1.0e-9 //1.0e-8 seems ok
#define TOL2 xxTOL //10 //1.e-8 // inextensibility constraint (more lax?)

#define broyden_maxiter 300

// Output (plotting/saving frequency)
#define plot_steps xxplot_steps
#define SaveExtendedSwimmerData xxSaveExtendedSwimmerData
#define SaveExtendedForcesData xxSaveExtendedForcesData

// Additional integer identifier to specify the simulation run
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define SimulationName "xxSimName" // "FreeNetwork-Nnetfil" STR(RandomMessNumNetworkFilaments) "-Kap" STR(Kap_NetworkFilaments)
#define SimulationTag "0000000"

#define PI  3.14159265358979323846264338327950288
#define PI2 6.28318530717959
#define TwoPI (2.*PI)

// Do not modify below here
/*****************************************************************************/


// Convenience (allow for multiple spellings) (as these are macros, this is OK)
#define OutputFolder                 "output/"
#define SimulationConfigName         OutputFolder SimulationName ".par"
#define SimulationTwoParticleConfigName       OutputFolder SimulationName "-2p.par"
#define SimulationDataName           OutputFolder SimulationName ".dat"
#define SimulationForceDataName           OutputFolder SimulationName "-forces.dat"
#define SimulationTwoParticleName    OutputFolder SimulationName "-2p.dat"
#define SimulationSwimmerDataName    OutputFolder SimulationName "-swimmer.dat"
#define SimulationGaitResetName      OutputFolder SimulationName "-gaitreset.dat"
#define SimulationGaitResetAndRemoveName      OutputFolder SimulationName "-gaitresetandremove.dat"
#define SimulationSwimmerDataCollisionName    OutputFolder SimulationName "-swimmer-collision.dat"
#define SimulationGaitResetCollisionName      OutputFolder SimulationName "-gaitreset-collision.dat"
#define SimulationSpringLinkDataName OutputFolder SimulationName "-springlinks.dat"
#define SimulationBackupDataName     OutputFolder SimulationName ".bak"

#define N_sw Nsw
#define N_worm Nworm
#define N_w Nworm
#define a_h ah
#define KAP Kap
#define SP4 Sp4
#define dl DL
#define dL DL
#define N_p Np
#define N_lam Nlam
#define N_broy Nbroy
#define mu MU
#define eta MU
#define timesteps TimeSteps
#define N Np



#endif
