#ifndef PROFILERS_INCLUDED
#define PROFILERS_INCLUDED

//#include "multi_filament_header.hpp"

struct timer {
	string name = "";
	string short_name = "";
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
	int value = 0;
	int count = 0; // times called

	void start();
	int end(); // Return time elapsed in milliseconds
	float average();
	void reset();
	void print_total();
	string print_total_short();
	void print_average();
};

struct profiler {

	timer totalTimeFCM = 			{.name = "this timestep: FCM           ",
									 .short_name = "  FCM   "};
	// (incl. collision barrier & 'omp parallel for' which can be slow.)
	timer totalTimeCopy = 			{.name = "this timestep: apply/assign F",
									 .short_name = "Assign F"};
	timer totalTimeBuildJacobian = 	{.name = "this timestep: build Jinv    ",
									 .short_name = "Bld Jinv"};
	timer totalTimeErrorCheck = 	{.name = "this timestep: error check   ",
									 .short_name = "ErrorChk"};
	timer totalTimeSolveJacobian = 	{.name = "this timestep: solve Jinv    ",
									 .short_name = "Slv Jinv"};
	timer totalTimestepTime = 		{.name = "this timestep: total         ",
									 .short_name = " Total  "};
	timer totalTime = 				{.name = "",
									 .short_name = "Time Left"};
	timer totalTimeCollisionBarrier={.name = "this timestep: collision bar ",
									 .short_name = "ColBarri"};

	timer fcm_make_zero = 			{.name = "fcm_make_zero                ",
									 .short_name = "FCMMakeZ"};
	timer fcm_gaussian_setup =		{.name = "fcm_gaussian_setup           ",
									 .short_name = "FCMGauss"};
	timer fcm_force_distribution =  {.name = "fcm_force_distribution       ",
								 	 .short_name = "FCMForcD"};
	timer fcm_fft_forward =  		{.name = "fcm_fft_forward              ",
 								 	 .short_name = "FCMFFTFw"};
	timer fcm_fft_backward =  		{.name = "fcm_fft_backward             ",
  								 	 .short_name = "FCMFFTBk"};
	timer fcm_flow_solve =  		{.name = "fcm_flow_solve               ",
   								 	 .short_name = "FCMFlowS"};
	timer fcm_particle_velocities_rotations = {.name = "fcm_particle_velocities_rotat",
								 	 .short_name = "FCMParVR"};

	void print_names();

 	void print_all(int nt);
};


#endif
