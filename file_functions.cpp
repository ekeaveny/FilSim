#include "file_functions.hpp"
#include "config.hpp" // Include parameters set by config.hpp

// Write simulation parameters to output file with ending ".par" (can add more if
// required)
void save_parameter_values_to_file(string parameter_file_name,
                                   double LfcmBox_x, double LfcmBox_y,
                                   double LfcmBox_z,
                                   ivec BeadNumbers, vec BendingFactors,
                                   vector<Filament>& filaments){

	// Only called by rank 0.

	ofstream OutputConfigFile(parameter_file_name, std::ios_base::out); //overwrites file
	OutputConfigFile<<"SimulationName"<<" "<<"SimulationTag"<<" "
	                <<"Nsw"<<" "<<"Nworm"<<" "<<"a"<<" "<<"ah"<<" "
	                <<"DL"<<" "<<"Sp4"<<" "<<"Kap"<<" "<<"C"<<" "<<"MU"<<" "
	                <<"curvature"<<" "<<"K0"<<" "<<"TimeSteps"<<" "
	                <<"StepsPerPeriod"<<" "<<"TOL"<<" "<<"gmres_tol"<<" "
	                <<"broyden_maxiter"<<" "<<"gmres_maxiter"<<" "
	                <<"LfcmBox_x"<<" "
	                <<"plot_steps"<<" "<<"STRAINTWISTX"<<" "
	                <<"STRAINTWISTY"<<" "<<"STRAINTWISTZ"<<" "
	                <<"Kap_NetworkFilaments"<<" "
	                <<"LfcmBox_z"<<" "<<"LfcmBox_y"<<" "<<"B"<<" "
	                <<"StepsPerSettlingTime"<<" "
	                <<"dt"<<" "
	                <<"initialFrictionSteps" << " "
	                <<"EnableSpringLinks" << " "
	                <<"ConnectedNodesNumSwimmingFilaments" << " "
	                <<"ConnectedNodesNumNodes" << " "
	                <<"ConnectedNodesConnectionRadius" << " "
	                <<"ConnectedNodesMinNodeSpacing" << " "
	                <<"ConnectedNodesSpacingAwayFromNode" << " "
	                <<"NetworkSpringConstant" << " "
	                <<"NetworkSpringNaturalLength" << " "
	                <<"\n";
	OutputConfigFile<<SimulationName<<" "<<SimulationTag<<" "
	                <<filaments.size()<<" "<<Nworm<<" "<<a<<" "<<ah<<" "
	                <<DL<<" "<<Sp4<<" "<<Kap<<" "<<C<<" "<<MU<<" "
	                <<curvature<<" "<<K0<<" "<<TimeSteps<<" "
	                <<StepsPerPeriod<<" "<<TOL<<" "<<0<<" "
	                <<broyden_maxiter<<" "<<0<<" "
	                <<LfcmBox_x<<" "
	                <<plot_steps<<" "<<STRAINTWISTX<<" "
	                <<STRAINTWISTY<<" "<<STRAINTWISTZ<<" "
	                <<Kap_NetworkFilaments<<" "
	                <<LfcmBox_z<<" "<<LfcmBox_y<<" "<<Bnumber<<" "
	                <<StepsPerSettlingTime<<" "
	                <<dt<<" "
	                <<initialFrictionSteps << " "
	                <<EnableSpringLinks << " "
	                <<ConnectedNodesNumSwimmingFilaments << " "
	                <<ConnectedNodesNumNodes << " "
	                <<ConnectedNodesConnectionRadius << " "
	                <<ConnectedNodesMinNodeSpacing << " "
	                <<ConnectedNodesSpacingAwayFromNode << " "
	                <<NetworkSpringConstant << " "
	                <<NetworkSpringNaturalLength << " "
	                <<endl<< endl;
	OutputConfigFile<<"BeadNumbers"<<endl;
	for (int i = 0; i < filaments.size(); i++) {
		OutputConfigFile << BeadNumbers[i];
		if (i < (filaments.size() - 1)) {
			OutputConfigFile << " ";
		}
	}
	OutputConfigFile << endl;
	OutputConfigFile << endl;

	OutputConfigFile<<"BendingModuli are KAP*"<<endl;
	for (int i = 0; i < filaments.size(); i++) {
		OutputConfigFile << BendingFactors[i];
		if (i < (filaments.size() - 1)) {
			OutputConfigFile << " ";
		}
	}
	OutputConfigFile << endl;

	OutputConfigFile.close();
}


// Write names of columns to simulation output/data file ".dat"
void save_data_column_names_to_file(string data_file_name,
                                    vector<Filament>& filaments){
	// Only called by rank 0.

	ofstream OutputFile(data_file_name,std::ios_base::out);
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	OutputFile << "it" << " ";
	for(int nn = 0; nn<filaments.size(); ++nn) { // do not attempt to parallelise this!
		filaments[nn].writeDataNames(OutputFile);
	}
	OutputFile << endl;
	OutputFile.close();
}


// Write names of columns to simulation output/data file ".dat"
void save_data_column_names_to_swimmer_velocity_file(string data_file_name,
                                                     vector<Filament>& filaments){
	// Only called by rank 0.

	ofstream OutputFile(data_file_name,std::ios_base::out);
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	OutputFile << "it" << " ";
	filaments[0].writeDataNames(OutputFile);
	filaments[0].writeVelocityDataNames(OutputFile);
	OutputFile << endl;
	OutputFile.close();
}

// Write names of columns to simulation output/data file ".dat"
void save_data_column_names_to_swimmer_collision_file(string data_file_name,
                                                      vector<Filament>& filaments){
	// Only called by rank 0.

	ofstream OutputFile(data_file_name,std::ios_base::out);
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	OutputFile << "it" << " ";
	filaments[0].writeCollisionDataNames(OutputFile);
	OutputFile << endl;
	OutputFile.close();
}

// Write names of columns to simulation output/data file ".dat"
void save_data_column_names_to_forces_file(string data_file_name,
                                           vector<Filament>& filaments){
	// Only called by rank 0.

	ofstream OutputFile(data_file_name,std::ios_base::out);
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	OutputFile << "it" << " ";
	filaments[0].writeForceDataNames(OutputFile);
	OutputFile << endl;
	OutputFile.close();
}


// Write position data to simulation output/data file ".dat"
void save_data_to_file(string data_file_name, vector<Filament>& filaments,
                       int nt){
	// Only called by rank 0.

	ofstream OutputFile (data_file_name,ios::app);
	OutputFile << nt << " ";
	for(int nn = 0; nn<filaments.size(); ++nn) { // do not attempt to parallelise!
		filaments[nn].writeData(OutputFile);
	}
	OutputFile << endl;
	OutputFile.close();
}

// Write position data to simulation output/data file ".dat"
void save_swimmer_velocity_data_to_file(string data_file_name,
                                        vector<Filament>& filaments,
                                        double **V,double **W, double **F,
                                        int nt){
	// Only called by rank 0.

	ofstream OutputFile (data_file_name,ios::app);
	OutputFile << nt << " ";
	filaments[0].writeData(OutputFile);
	filaments[0].writeVelocityData(OutputFile, V, W, F, nt);
	OutputFile << endl;
	OutputFile.close();
}

// Write force data to simulation output/data file ".dat"
void save_force_data_to_file(string data_file_name,
                             vector<Filament>& filaments,
                             double **F, double **T,
                             int nt){
	// Only called by rank 0.

	ofstream OutputFile (data_file_name,ios::app);
	OutputFile << nt << " ";
	for(int nn = 0; nn<filaments.size(); ++nn) { // do not attempt to parallelise!
		filaments[nn].writeForceData(OutputFile, F, T, nt);
	}
	OutputFile << endl;
	OutputFile.close();
}

// Write collision data to simulation output/data file ".dat"
void save_swimmer_collision_data_to_file(string data_file_name,
                                         vector<Filament>& filaments,
                                         double F[][3],
                                         int nt){
	// Only called by rank 0.

	ofstream OutputFile (data_file_name,ios::app);
	OutputFile << nt << " ";
	filaments[0].writeCollisionData(OutputFile, F, nt);
	OutputFile << endl;
	OutputFile.close();
}

// Write position/orientation data to simulation backup file ".bak" which
// contains X, U, q data from the last two timesteps. This way if the sim
// needs restarting, the data for the timestepping/initial guessing is there.
void save_backup_data_to_file(string backup_data_file_name,
                              vector<Filament>& filaments,
                              int nt){
	// Only called by rank 0.

	ofstream OutputFile (backup_data_file_name); // Overwrite, don't append.
	OutputFile << nt << " ";
	for(int nn = 0; nn<filaments.size(); ++nn) { // do not attempt to parallelise!
		filaments[nn].writeBackupData(OutputFile);
	}
	OutputFile.close();
}

// Load nt, X, q, U, Xt, qt, Ut, lam, lam1, lam2 from file and
// into  filaments .
// NOTE: Not set up for filaments of different lengths or beads of different
// sizes.
tuple<int, int, int> continue_from_checkpoint(string data_file_name,
                                              vector<Filament>& filaments){
	// Only called by rank 0.

	string backup_data_file_name = OutputFolder + data_file_name + ".bak";
	string parameter_data_file_name = OutputFolder + data_file_name + ".par";
	int nt = -1;
	int forceArrayLengthFilled = 0;
	int stateArrayLengthFilled = 0;
	int N_sw_override;

	// Input number of filaments and filament lengths
	ifstream input_file (parameter_data_file_name);
	try {
		string s, num_fils_as_string;
		getline(input_file,s); // First \n-delimited item. Ignore
		getline(input_file,s,' '); // 1st space delimited item. Ignore
		getline(input_file,s,' '); // 2nd space delimited item. Ignore
		getline(input_file,num_fils_as_string,' '); // 3rd space delimited item.
		cout << "N_sw_override is..." << endl;
		N_sw_override = stoi(num_fils_as_string);
		cout << "..." << N_sw_override << endl;

		getline(input_file,s); // Rest of 2nd line. Ignore
		getline(input_file,s); // 3rd line. Ignore
		getline(input_file,s); // 4th line. Ignore

		// bool legacy_vertical_num_beads = false;
		for (int i=0; i<N_sw_override; ++i) {
			string num_beads_string;
			// If STOI error, probably because of this line; probably because
			// list of beads is vertical (\n separated) rather than horizontal
			// (space separated), like in the newer .par files.
			// if (legacy_vertical_num_beads == false){
			//     getline(input_file,num_beads_string,' ');
			// } else {
			//     getline(input_file,num_beads_string);
			// }
			// if (num_beads_string == ""){
			//     legacy_vertical_num_beads = true;
			//     getline(input_file,num_beads_string);
			// }
			getline(input_file,num_beads_string,' ');

			// cout << "num_beads is..." << endl;
			// cout << num_beads_string << " which is " << endl;
			// cout << stoi(num_beads_string) << "." << endl;

			int num_beads = stoi(num_beads_string);
			Filament new_filament(num_beads);
			filaments.push_back(new_filament);
		}
	}
	catch (const std::exception& e) {
		cout << endl << "ERROR: Parameter file '" << parameter_data_file_name
		     << "' either does not exist or is corrupted." << endl << endl;
		cout << "Exception raised: " << e.what() << endl;
		exit(-1);
	}
	input_file.close();


	// Input positions from backup file
	ifstream input_file2 (backup_data_file_name);
	try {
		string nt_as_string;
		getline(input_file2,nt_as_string,' '); // First space-delimited item
		nt = stoi(nt_as_string);
		for (int i=0; i<N_sw_override; ++i) {
			filaments[i].continue_from_checkpoint(input_file2,
			                                      forceArrayLengthFilled,
			                                      stateArrayLengthFilled);

			forceArrayLengthFilled += 1*filaments[i].length();
			stateArrayLengthFilled += 6*filaments[i].length();
		}
	}
	catch (const std::exception& e) {
		cout << endl << "ERROR: Checkpoint file '" << backup_data_file_name
		     << "' either does not exist or is corrupted." << endl << endl;
		exit(-1);
	}
	input_file2.close();

	// compute global number of particles
	int Np = 0;
	for(int nn = 0; nn<filaments.size(); ++nn) {
		Np += filaments[nn].length();
	}

	return make_tuple(nt+1, forceArrayLengthFilled, stateArrayLengthFilled);
}
