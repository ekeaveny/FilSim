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
#define SwimProblem true
#define SedimentationProblem false

// Initial setup
#define GTinitialisation false
#define RANDinitialisation false
#define StackedInitialisation false
#define GTTallBoxInitialisation false

#define RandomMessToSwimThroughInitialisation true
#define RandomMessNumNetworkFilaments 10
#define RandomMessNumSwimmingFilaments 1
// Will put centre of filament inside a box of size
// [0,RandomMessBoxSize_xy]*[0,RandomMessBoxSize_xy]*[0,RandomMessBoxSize_z]
#define RandomMessBoxSize_x 77.6 // 38.8 // 116.412000258 // 465.648 // 38.8 // 77.60
#define RandomMessBoxSize_y 38.8 // 116.412000258 // 465.648 // 38.8 // 77.60
                               // 77.608 equivalent to box with NPTS=256. 116.412 equivalent to NPTS=384
#define RandomMessBoxSize_z 38.8 // 465.648001033


// Number of filaments/swimmers
// in the present implementation, this number is limited by the stack size, if ~STACKALLOCATION
// For lattices, note that for Nx,Ny,Nz grid points, num of filaments required = 3 Nx Ny Nz - Nx Ny - Ny Nz - Nx Nz.
//#define Nsw (3*LatticeNi*LatticeNj*LatticeNk - LatticeNi*LatticeNj - LatticeNj*LatticeNk - LatticeNj*LatticeNi + LatticeNumFilamentsOutside) // = 54 (for lattice) + 1
#define Nsw 11 // 401 // (RandomMessNumNetworkFilaments+RandomMessNumSwimmingFilaments)

// Number of particles per filament/swimmer
#define Nworm 15
#define max_Nworm 15 // Maximum possible filament length (for array allocation)

// Number of filaments/swimmers
//#define Nsw 10
// Number of particles per filament/swimmer
//#define Nworm 20

// FCM configuration
#define NPTS_X 256L // 384L // 128 // 256
#define NPTS_Y 128L // 384L // 128 // 256
#define NPTS_Z 128L//1536L // 128 // 256
#define ngd 20//15

// constants
#define a  1.0 // segment radius
#define ah 1.0 // head radius
#define DL 2.2 // segment distance
#define MU 1.0 // viscsoty
#define WeightPerLengthX 0.0 // gravitational force line density in +X-direction
#define WeightPerLengthY 0.0 // gravitational force line density in +Y-direction
#define WeightPerLengthZ -1.0 // gravitational force line density in +Z-direction
// #define Lbox 2. // initial box size in swimmer lengths
#define InitialSeparation 50. //separation in GT initialisation

// Sperm number and bending modulus
#define Sp4       24.0 // not Sp, but Sp^4
#define Bnumber   1  //1.e-3 // elasto-gravitational number

#if SwimProblem
    #define Kap                   1000    //1800.0 // in terms of Bnumber, Kap = WL^3/B
    #define C                     Kap     //1800.0
    #define Kap_NetworkFilaments  3.6      // (Kap=36 corresponds to B=1000. Kap=36000 corresponds to B=1. Of course B is not really well defined without gravity.)
    #define C_NetworkFilaments    Kap_NetworkFilaments
#endif

// Target curvature per length (double!)
#define curvature 8.25

// Swimming: Wave number, helical swimming alpha and beta
#define K0 1.
#define SwimmingHelixAlpha 0.
#define SwimmingHelixBeta 1.

// Time integration
#define StepsPerPeriod 3000  // for swimmer problem
                            // note StepsPerPeriod should be a multiple of 100 unless you redefine plot_steps below.
#define StepsPerSettlingTime 0 // for sedimentation problem
#define TimeSteps 180000//      StepsPerPeriod/30 // 30000

// Strain twist vectors
#define STRAINTWISTX 0.
#define STRAINTWISTY 0.
#define STRAINTWISTZ 0.

// number of initial time steps computed with implicit euler instead of BDF scheme
#define ImplEulerSteps 1//2//100000 //1
#define initialFrictionSteps 1 // number of steps with local friction only

// Tolerance for non-linear system, gmres, and maximum number of iterations
#define TOL 1.e-8//10 //1.e-8 //1.0e-9 //1.0e-8 seems ok
#define TOL2 1.e-8//10 //1.e-8 // inextensibility constraint (more lax?)

#define broyden_maxiter 100

// Output (plotting/saving frequency)
#define plot_steps 3

// Additional integer identifier to specify the simulation run
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define SimulationName "Test247" // "Swim-Nnetfil400-Kap3p6-NPTSxy384z1536" // "FreeNetwork-Nnetfil" STR(RandomMessNumNetworkFilaments) "-Kap" STR(Kap_NetworkFilaments)
#define SimulationTag "0000000"

#define PI  3.14159265358979323846264338327950288
#define PI2 6.28318530717959
#define TwoPI (2.*PI)

// Do not modify below here
/*****************************************************************************/


// Convenience (allow for multiple spellings) (as these are macros, this is OK)
#define OutputFolder               "output/"
#define SimulationConfigName       OutputFolder SimulationName ".par"
#define SimulationDataName         OutputFolder SimulationName ".dat"
#define SimulationSwimmerDataName  OutputFolder SimulationName "-swimmer.dat"
#define SimulationBackupDataName   OutputFolder SimulationName ".bak"

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
