#include "print_functions.hpp"
#include "config.hpp" // Include parameters set by config.hpp

// Formats a time in seconds as something more readable and
// (mostly) fixed width.
//
//       format_time(20)     returns '   20.0s'
//       format_time(200)    returns '    3:20'
//       format_time(20000)  returns ' 5:33:20'
//       format_time(200000) returns '2d 5:33:20'
string format_time(float elapsed_time){
	int tr2b_d, tr2b_h, tr2b_m, tr2b_s;
	char* elapsed_time_hms = new char[20];
	if (elapsed_time >= 86400) {
		tr2b_m = floor(elapsed_time/60);
		tr2b_s = fmod(elapsed_time,60);
		tr2b_h = floor(tr2b_m/60);
		tr2b_m = fmod(tr2b_m,60);
		tr2b_d = floor(tr2b_h/24);
		tr2b_h = fmod(tr2b_h,24);
		sprintf(elapsed_time_hms,"%dd%2d:%02d:%02d", tr2b_d, tr2b_h, tr2b_m, tr2b_s);
	}
	else if (elapsed_time >= 3600) {
		tr2b_m = floor(elapsed_time/60);
		tr2b_s = fmod(elapsed_time,60);
		tr2b_h = floor(tr2b_m/60);
		tr2b_m = fmod(tr2b_m,60);
		sprintf(elapsed_time_hms,"%2d:%02d:%02d", tr2b_h, tr2b_m, tr2b_s);
	}
	else if (elapsed_time >= 60) {
		tr2b_m = floor(elapsed_time/60);
		tr2b_s = fmod(elapsed_time,60);
		sprintf(elapsed_time_hms,"   %2d:%02d",tr2b_m, tr2b_s);
	}
	else {
		sprintf(elapsed_time_hms,"%7.1fs",elapsed_time);
	}
	std::string str = elapsed_time_hms;
	delete[] elapsed_time_hms;
	return str;
}


// Prints the simulation parameter values in a pretty table.
void print_parameter_values_in_table(int Np, double LfcmBox_x, double LfcmBox_y,
                                     double LfcmBox_z, double L,
                                     double omega, double dt){

	double volfrac = 4.0*a*a*a*PI/3.0;
	volfrac *= (N)/(LfcmBox_x*LfcmBox_y*LfcmBox_z)*100;

	cout << "----------------------------------------------------------------" << endl;
	cout << "                     SIMULATION PARAMETERS" << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << " Np                                 " << Np              << endl;
	cout << " N_sw                               " << N_sw            << endl;
	cout << " N_w                                " << N_w             << endl << endl;

	cout << " Fluid solver                       FCM (periodic)"     << endl;
	cout << " NPTS_X                             " << NPTS_X         << endl;
	cout << " NPTS_Y                             " << NPTS_Y         << endl;
	cout << " NPTS_Z                             " << NPTS_Z         << endl;
	cout << " LfcmBox_x                          " << LfcmBox_x      << endl;
	cout << " LfcmBox_y                          " << LfcmBox_y      << endl;
	cout << " LfcmBox_z                          " << LfcmBox_z      << endl;
	cout << " (volume fraction)                  " << volfrac <<"%"   << endl << endl;

	cout << " L (NB fil length might be changed) " << L               << endl;
	cout << " DL                                 " << DL              << endl;

	cout << "----------------------------------------------------------------" << endl;
	#if SedimentationProblem
		cout << " B                                  " << Bnumber         << endl;
		cout << " Weight/L                           " << WeightPerLengthX << ", " << WeightPerLengthY << ", " << WeightPerLengthZ << endl;
		cout << " Kap                                " << Kap             << endl << endl;
	#endif
	#if RandomMessToSwimThroughInitialisation || RandomMessToSwimThroughInitialisationAllowOverlap
		cout << " RandomMessToSwimThroughInitialisation specifics:" << endl;
		cout << endl;
		cout << " NETWORK FILAMENTS" << endl;
		cout << "   RandomMessNumNetworkFilaments    " << RandomMessNumNetworkFilaments << endl;
		cout << "   Centres of filaments inside box  [0," << RandomMessBoxSize_x << "]*[0," << RandomMessBoxSize_y << "]*[0," << RandomMessBoxSize_z << "]" << endl;
		cout << "   Kap_NetworkFilaments             " << Kap_NetworkFilaments << endl;
		cout << "   C_NetworkFilaments               " << C_NetworkFilaments << endl;
		cout << " SWIMMERS" << endl;
		cout << "   RandomMessNumSwimmingFilaments   " << RandomMessNumSwimmingFilaments << endl;
		cout << "   Kap                              " << Kap << endl;
		cout << "   C                                " << C << endl;
		cout << "   omega = 2πf = Kap*Sp^4/(mu*L^4)  " << omega << endl;
		cout << "   Sp^4                             " << Sp4 << endl;
	#endif
	cout << "----------------------------------------------------------------" << endl;
	cout << " Collision barrier FS over 2a       " << Barrier_FS_over_2a       << endl;
	cout << " Collision barrier distance cap     " << Barrier_distance_cap     << endl;	
	cout << "----------------------------------------------------------------" << endl;
	cout << " timesteps                          " << timesteps       << endl;
	cout << " dt = 1/(StepsPerPeriod*f)          " << dt << endl;
	cout << " (final time)                       " << dt*timesteps << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << " Simulation Name                    " << SimulationName << endl;
	cout << " Simulation Tag                     " << SimulationTag << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << endl;
}


// Sets the number of significant figures for screen output
void set_precision_of_screen_output(int precision){
	cout.precision(precision);
	cout << std::scientific;
}
