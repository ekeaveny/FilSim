#include "filament_initialisation_functions.hpp"
#include "spring_link_functions.hpp"
#include "config.hpp" // Include parameters set by config.hpp

/**
    filament_position_initialisation   Initialises positions of filaments, as
                                       given by the switches in config.hpp.
                                       NOTE: This function alters the vector
                                       `filaments` .

                                       Properties of the filaments are imposed
                                       by `filament_property_initialisation`
 */
tuple<int, int> filament_position_initialisation(vector<Filament>& filaments
#if EnableSpringLinks
	                                                 , spring_links& spring_links
#endif
                                                 ) {

	vec HeadTangent(3,fill::zeros);  HeadTangent(0)  = 1.;
	vec HeadNormal(3,fill::zeros);   HeadNormal(1)   = 1.;
	vec HeadBiNormal(3,fill::zeros); HeadBiNormal(2) = 1.;
	vec HeadPosition(3);

	// so that filament knows where it sits in the force vectors etc.
	int forceArrayLengthFilled = 0;
	// so that filament knows where it sits in the state vectors etc.
	int stateArrayLengthFilled = 0;

	// -------------------------------------------------------------------------

	#if GTinitialisation
		// Initialisation in  a cylinder, see Gustavsson and Tornberg paper.

		for(int nn = 0; nn<Nsw; ++nn) {
			bool Overlap(true);
			do {
				// Arbitrary stack of swimmers in normals dir.
				// (Replace this for random positions and orientations)
				HeadPosition = .5 * InitialSeparation
				               * (  std::cos(2*PI*nn/Nsw)*HeadNormal
				                    + std::sin(2*PI*nn/Nsw)*HeadBiNormal );
				cout << "Proposing head position ... \n" << endl
				     << HeadPosition << endl << endl;
				filaments[nn].initialFilamentSetup(
					HeadPosition, HeadTangent, HeadNormal,
					forceArrayLengthFilled, stateArrayLengthFilled);
				Overlap = filaments[nn].checkOverlapWith(filaments,nn);
			}
			while(Overlap);
			cout << "... accepted head position."<<endl<<endl;
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		// ---------------------------------------------------------------------

	#elif RANDinitialisation
		// random cloud of particles

		int randSeed = 2;
		arma_rng::set_seed(randSeed);

		for(int nn = 0; nn<Nsw; ++nn) {
			bool Overlap(true);
			do {
				cout << "HeadTangent "<< HeadTangent << endl;
				cout << "HeadNormal "<< HeadNormal << endl;
				HeadPosition.randu();
				HeadPosition(0) *= .9*LfcmBox_x;
				HeadPosition(1) *= .9*LfcmBox_y;
				HeadPosition(2) *= .9*LfcmBox_z;
				cout << "Proposing head position ... \n" << endl
				     << HeadPosition << endl << endl;
				filaments[nn].initialFilamentSetupIsoRand(
					HeadPosition,
					forceArrayLengthFilled, stateArrayLengthFilled);
				Overlap = filaments[nn].checkOverlapWith(filaments,nn);
			}
			while(Overlap);
			cout << "... accepted head position."<<endl<<endl;
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		// ---------------------------------------------------------------------

	#elif RandomMessToSwimThroughInitialisation
		// Random mess to swim through. Seed randomly and don't let anything
		// overlap.

		float x_coord; float y_coord; float z_coord;
		int nn = 0;
		float kap_in; float c_in;

		// Put swimmer in first, at top, tail in the air.
		for(int nn = 0; nn<RandomMessNumSwimmingFilaments; ++nn) {
			Filament new_filament(Nworm_swimmer);
			filaments.push_back(new_filament);
			x_coord = LfcmBox_x/2; // 0.1*10;
			y_coord = LfcmBox_y/2;
			z_coord = 0; // LfcmBox_z/2;
			HeadPosition = {x_coord, y_coord, z_coord};
			HeadTangent = {0.,0.,1. }; // Pointing towards tail
			HeadNormal = {1.,0.,0. };
			filaments[nn].initialFilamentSetup(
				HeadPosition, HeadTangent, HeadNormal,
				forceArrayLengthFilled,  stateArrayLengthFilled);
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << RandomMessNumSwimmingFilaments
		     << " swimming worms placed." << endl;

		// Network of filaments
		//int randSeed = 3;
		// randSeed now set in config
		arma_rng::set_seed(randSeed);
		int nn_local = 0;
		for(nn = RandomMessNumSwimmingFilaments; nn<N_sw; ++nn) {
			Filament new_filament(Nworm);
			filaments.push_back(new_filament);
			nn_local = nn - RandomMessNumSwimmingFilaments;
			bool Overlap(true);
			cout << "Proposing CoM position of network filament "
			     << nn_local+1 << "/" << RandomMessNumNetworkFilaments
			     << "...";
			do {
				// Uniform random in box of size RandomMessBoxSize
				HeadPosition.randu(); // Actually CoM position, not head pos'n
				HeadPosition[0] *= RandomMessBoxSize_x;
				HeadPosition[1] *= RandomMessBoxSize_y;
				HeadPosition[2] *= RandomMessBoxSize_z;
				filaments[nn].initialFilamentSetupIsoRandCentred(
					HeadPosition,
					forceArrayLengthFilled, stateArrayLengthFilled);
				Overlap = filaments[nn].checkOverlapWithPeriodic(filaments,
				                                                 nn, 1.05);
				//                                                                                                   was 2
			}
			while(Overlap);
			cout << " done." << endl;
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << RandomMessNumNetworkFilaments
		     << " filaments initialised randomly." << endl;
		cout << N_sw << " filaments initialised in total." <<endl;

		// ---------------------------------------------------------------------

	#elif RandomMessToSwimThroughInitialisationAllowOverlap
		// Random mess to swim through. Seed randomly and let them Overlap
		// (well, not directly but pretty close.) This is Part 1 of a two-
		// stage algorithm for creating dense networks. Seed like this, then
		// add all forces, but with Stokes drag only. Run until overlap > 2.

		float x_coord; float y_coord; float z_coord;
		int nn = 0;
		float kap_in; float c_in;

		// Put swimmer in first, at top, tail in the air.
		for(int nn = 0; nn<RandomMessNumSwimmingFilaments; ++nn) {
			Filament new_filament(Nworm_swimmer);
			filaments.push_back(new_filament);
			x_coord = LfcmBox_x/2; // 0.1*10;
			y_coord = LfcmBox_y/2;
			z_coord = 0; // LfcmBox_z/2;
			HeadPosition = {x_coord, y_coord, z_coord};
			HeadTangent = {0.,0.,1. }; // Pointing towards tail
			HeadNormal = {1.,0.,0. };
			filaments[nn].initialFilamentSetup(
				HeadPosition, HeadTangent, HeadNormal,
				forceArrayLengthFilled,  stateArrayLengthFilled);
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << RandomMessNumSwimmingFilaments
		     << " swimming worms placed." << endl;

		// Network of filaments
		//int randSeed = 3;
		// randSeed now set in config
		arma_rng::set_seed(randSeed);
		int nn_local = 0;
		for(nn = RandomMessNumSwimmingFilaments; nn<N_sw; ++nn) {
			Filament new_filament(Nworm);
			filaments.push_back(new_filament);
			nn_local = nn - RandomMessNumSwimmingFilaments;
			bool Overlap(true);
			cout << "Proposing CoM position of network filament "
			     << nn_local+1 << "/" << RandomMessNumNetworkFilaments
			     << "...";
			do {
				// Uniform random in box of size RandomMessBoxSize
				HeadPosition.randu(); // Actually CoM position, not head pos'n
				HeadPosition[0] *= RandomMessBoxSize_x;
				HeadPosition[1] *= RandomMessBoxSize_y;
				HeadPosition[2] *= RandomMessBoxSize_z;
				filaments[nn].initialFilamentSetupIsoRandCentred(
					HeadPosition,
					forceArrayLengthFilled, stateArrayLengthFilled);
				Overlap = filaments[nn].checkOverlapWithPeriodic(filaments,
				                                                 nn, 0.01);
			}
			while(Overlap);
			cout << " done." << endl;
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << RandomMessNumNetworkFilaments
		     << " filaments initialised randomly, allowed to overlap." << endl;
		cout << N_sw << " filaments initialised in total, allowed to overlap." <<endl;

		// ---------------------------------------------------------------------

	#elif RandomMessToSwimThroughInitialisationTwoParticle
		// Random mess to swim through. Seed randomly and don't let anything
		// overlap.

		float x_coord; float y_coord; float z_coord;
		int nn = 0;
		float kap_in; float c_in;

		// Network of filaments
		//int randSeed = 3;

		int nn_local = 0;
		for(nn = 0; nn<N_sw-2; ++nn) {
			Filament new_filament(Nworm);
			filaments.push_back(new_filament);
			nn_local = nn;// - RandomMessNumSwimmingFilaments;
			bool Overlap(true);
			cout << "Proposing CoM position of network filament "
			     << nn_local+1 << "/" << RandomMessNumNetworkFilaments
			     << "...";
			do {
				// Uniform random in box of size RandomMessBoxSize
				HeadPosition.randu(); // Actually CoM position, not head pos'n
				HeadPosition[0] *= RandomMessBoxSize_x;
				HeadPosition[1] *= RandomMessBoxSize_y;
				HeadPosition[2] *= RandomMessBoxSize_z;
				filaments[nn].initialFilamentSetupIsoRandCentred(
					HeadPosition,
					forceArrayLengthFilled, stateArrayLengthFilled);
				Overlap = filaments[nn].checkOverlapWithPeriodic(filaments,
				                                                 nn, 1.05);
				//                                                                                                   was 2
			}
			while(Overlap);
			cout << " done." << endl;
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << RandomMessNumNetworkFilaments
		     << " filaments initialised randomly." << endl;
		cout << N_sw-2 << " filaments initialised in total." <<endl;

		// ---------------------------------------------------------------------

	#elif ConnectedNodesInitialisation || ConnectedNodesInitialisationAllowOverlap
		// Randomly seed nodes. Then connect them with filaments.
		// Filaments are therefore of different lengths.
		// These filaments can be placed, or not placed, based on a probability.

		float x_coord, y_coord, z_coord;
		float kap_in, c_in;
		vec NodePosition(3);
		vector<vec> Nodes;
		vector<vector<int> > node_to_node_filaments;
		int num_node_to_node_filaments_suggested, N_sw_override;
		int max_fil_length_so_far = 0;

		// Put swimmer in first, at top, tail in the air.
		for(int nn = 0; nn<ConnectedNodesNumSwimmingFilaments; ++nn) {
			max_fil_length_so_far = max(max_fil_length_so_far,
			                            Nworm_swimmer);
			Filament new_filament(Nworm_swimmer);
			filaments.push_back(new_filament);
			x_coord = LfcmBox_x/2; // 0.1*10;
			y_coord = LfcmBox_y/2;
			z_coord = 0; // LfcmBox_z/2;
			HeadPosition = {x_coord, y_coord, z_coord};
			HeadTangent = {0.,0.,1. }; // Pointing towards tail
			HeadNormal = {1.,0.,0. };
			filaments[nn].initialFilamentSetup(
				HeadPosition, HeadTangent, HeadNormal,
				forceArrayLengthFilled,  stateArrayLengthFilled);
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << ConnectedNodesNumSwimmingFilaments
		     << " swimming worms placed." << endl;

		// Put in nodes
		arma_rng::set_seed(randSeed); // randSeed set in config

		for(int nn = 0; nn < ConnectedNodesNumNodes; ++nn) {
			bool Overlap(true);
			cout << "Proposing node position "
			     << nn+1 << "/" << ConnectedNodesNumNodes
			     << "..." << endl;
			Nodes.push_back({0,0,0});
			do {
				// Uniform random in box of size RandomMessBoxSize
				NodePosition.randu();
				NodePosition[0] *= ConnectedNodesBoxSize_x;
				NodePosition[1] *= ConnectedNodesBoxSize_y;
				NodePosition[2] *= ConnectedNodesBoxSize_z;
				Nodes[nn] = NodePosition;
				Overlap = checkNodeOverlapWithPeriodic(Nodes, nn, ConnectedNodesMinNodeSpacing);
			}
			while(Overlap);
			//cout << " done." << endl;
			cout << NodePosition[0] << ", " << NodePosition[1] << ", " << NodePosition[2] << endl;
		}
		cout << ConnectedNodesNumNodes
		     << " nodes placed." << endl;

		//cin.get();

		spring_links.set_num_nodes(Nodes.size());

		// Which nodes should the network filaments go between?
		node_to_node_filaments = decide_node_to_node_filaments(Nodes);
		num_node_to_node_filaments_suggested = node_to_node_filaments.size();
		cout << num_node_to_node_filaments_suggested
		     << " network filaments suggested." << endl;

		if (num_node_to_node_filaments_suggested>Nsw) {
			cout << endl;
			cout << "ERROR: Nsw < NUMBER OF FILAMENTS SUGGESTED." << endl;
			cout << "SET Nsw TO MAXIMUM NUMBER OF FILAMENTS." << endl;
			//cout << "SEGFAULT IMMINENT." << endl;
			cout << endl;
			exit (EXIT_FAILURE);
		}

		// Place in the filaments. If it overlaps, don't place it, just forget it.
		//N_sw_override = ConnectedNodesNumSwimmingFilaments
		//                + num_node_to_node_filaments;
		int num_node_to_node_filaments_placed = 0;

		for(int i = 0; i<num_node_to_node_filaments_suggested; ++i) {
			int nn = (ConnectedNodesNumSwimmingFilaments
			          + num_node_to_node_filaments_placed);
			bool Overlap(true);
			cout << "Proposing position of suggested network filament "
			     << i+1 << "/" << num_node_to_node_filaments_suggested
			     << "...";

			vec node1 = Nodes[node_to_node_filaments[i][0]];
			vec node2 = Nodes[node_to_node_filaments[i][1]];

			int num_beads_required;
			tie(num_beads_required, ignore, ignore) =
				beads_between_nodes(node1, node2);
			Filament new_filament(num_beads_required);
			//cout << "X" << num_beads_required << endl;
			new_filament.initialFilamentSetupBetweenTwoPoints(
				node1,
				node2,
				forceArrayLengthFilled, stateArrayLengthFilled);
			//cout << "Y" << endl;
			#if ConnectedNodesInitialisationAllowOverlap
				Overlap = new_filament.checkOverlapWithPeriodic(filaments,
				                                                nn, 0.01);
			#else
				Overlap = new_filament.checkOverlapWithPeriodic(filaments,
				                                                nn, 1.05);
			#endif
			//cout << "Z" << endl;


			if (!Overlap) {
				max_fil_length_so_far = max(max_fil_length_so_far,
				                            filaments[nn].length());
				filaments.push_back(new_filament);
				cout << " done (length " << filaments[nn].length() << ") "
				     << "(becoming global fil. no. " << nn+1 << ")" << endl;
				forceArrayLengthFilled += 1*filaments[nn].length();
				stateArrayLengthFilled += 6*filaments[nn].length();
				++num_node_to_node_filaments_placed;

				// Add {fil,bead} to filaments_at_nodes[node]
				spring_links.filaments_at_nodes[node_to_node_filaments[i][0]].push_back({nn,0});
				spring_links.filaments_at_nodes[node_to_node_filaments[i][1]].push_back(
					{nn,filaments[nn].length()-1}
					);
			} else {
				cout << " failed, due to overlap. Skipping." << endl;
			}
		}

		#if ConnectedNodesInitialisationAllowOverlap
			cout << num_node_to_node_filaments_placed
			     << " network filaments initialised, allowing overlap." << endl;
		#else
			cout << num_node_to_node_filaments_placed
			     << " network filaments initialised." << endl;
		#endif

		N_sw_override = (ConnectedNodesNumSwimmingFilaments
		                 + num_node_to_node_filaments_placed);
		cout << N_sw_override << " filaments initialised in total." <<endl;
		cout << "Maximum filament length: " << max_fil_length_so_far
		     << " segments." << endl;


		// The springs between the filaments at each node are created elsewhere.

		// ---------------------------------------------------------------------

	#elif GTTallBoxInitialisation
		// Tall box for sedimenting filaments in periodic domain.
		// Inspired by G&T fig 5.

		float x_coord; float y_coord; float z_coord;
		int nn = 0;
		float kap_in; float c_in;

		// Network of filaments
		int randSeed = 3;
		arma_rng::set_seed(randSeed);
		int nn_local = 0;
		for(nn = 0; nn<N_sw; ++nn) {
			nn_local = nn;
			bool Overlap(true);
			cout << "Proposing CoM position of filament "
			     << nn_local+1 << "/" << N_sw
			     << "...";
			do {
				// Uniform random in box of size RandomMessBoxSize
				HeadPosition.randu(); // Actually CoM position, not head pos'n
				HeadPosition[0] *= LfcmBox_x;
				HeadPosition[1] *= LfcmBox_y;
				HeadPosition[2] *= LfcmBox_z;
				filaments[nn].initialFilamentSetupIsoRandCentred(
					HeadPosition,
					forceArrayLengthFilled, stateArrayLengthFilled);
				Overlap = filaments[nn].checkOverlapWithPeriodic(filaments,
				                                                 nn, 2);
			}
			while(Overlap);
			cout << " done." << endl;
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << N_sw << " filaments initialised in total." <<endl;

		// ---------------------------------------------------------------------

	#elif TempInitialisation

		float x_coord; float y_coord; float z_coord;
		int nn = 0;
		float kap_in; float c_in;

		// Put swimmer in first
		for(int nn = 0; nn<RandomMessNumSwimmingFilaments; ++nn) {
			x_coord = LfcmBox_x/2; // 0.1*10;
			y_coord = LfcmBox_y/2;
			z_coord = LfcmBox_z/2;
			HeadPosition = {x_coord, y_coord, z_coord};
			HeadTangent = {-1.,0.,0. }; // Pointing towards tail
			HeadNormal = {0.,-1.,0. };
			filaments[nn].initialFilamentSetup(
				HeadPosition, HeadTangent, HeadNormal,
				forceArrayLengthFilled,  stateArrayLengthFilled);
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << RandomMessNumSwimmingFilaments
		     << " swimming worms placed." << endl;

		// Network of filaments
		int randSeed = 3;
		arma_rng::set_seed(randSeed);
		int nn_local = 0;
		for(nn = RandomMessNumSwimmingFilaments; nn<N_sw; ++nn) {
			nn_local = nn - RandomMessNumSwimmingFilaments;
			x_coord = LfcmBox_x/2; // 0.1*10;
			y_coord = LfcmBox_y/2;
			z_coord = LfcmBox_z/2 + 2.2;
			HeadPosition = {x_coord, y_coord, z_coord};
			HeadTangent = {-1.,0.,0. }; // Pointing towards tail
			HeadNormal = {0.,-1.,0. };
			filaments[nn].initialFilamentSetup(
				HeadPosition, HeadTangent, HeadNormal,
				forceArrayLengthFilled,  stateArrayLengthFilled);
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << RandomMessNumNetworkFilaments
		     << " filaments initialised randomly." << endl;
		cout << N_sw << " filaments initialised in total." <<endl;

	#elif StackedInitialisation

		float x_coord; float y_coord; float z_coord;
		float kap_in; float c_in;

		// Put swimmer in first
		for(int nn = 0; nn<N_sw; ++nn) {
			x_coord = LfcmBox_x/2; // 0.1*10;
			y_coord = LfcmBox_y/2;
			z_coord = LfcmBox_z/2;
			HeadPosition = {x_coord, y_coord, z_coord + 3*nn};
			HeadTangent = {-1.,0.,0. }; // Pointing towards tail
			HeadNormal = {0.,-1.,0. };
			filaments[nn].initialFilamentSetup(
				HeadPosition, HeadTangent, HeadNormal,
				forceArrayLengthFilled,  stateArrayLengthFilled);
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}
		cout << N_sw
		     << " worms placed." << endl;

	#elif TestInitialisation1

		float x_coord, y_coord, z_coord, LL;
		LL = 15*2.2;
		x_coord = LfcmBox_x/2 + 0.3*LL; // 0.1*10;
		y_coord = LfcmBox_y/2;
		z_coord = LfcmBox_z/2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {-1.,0.,0. }; // Pointing towards tail
		HeadNormal = {0.,-1.,0. };
		filaments[0].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[0].length();
		stateArrayLengthFilled += 6*filaments[0].length();

		x_coord = LfcmBox_x/2 + 0.7*LL; // 0.1*10;
		y_coord = LfcmBox_y/2 + 0.5*LL;
		z_coord = LfcmBox_z/2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {0.,-1.,0. }; // Pointing towards tail
		HeadNormal = {0.,0.,-1. };
		filaments[1].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[1].length();
		stateArrayLengthFilled += 6*filaments[1].length();

	#elif TestInitialisation2

		float x_coord, y_coord, z_coord, LL;
		LL = 15*2.2;
		x_coord = LfcmBox_x/2 - 0.7*LL; // 0.1*10;
		y_coord = LfcmBox_y/2;
		z_coord = LfcmBox_z/2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {1.,0.,0. }; // Pointing towards tail
		HeadNormal = {0.,1.,0. };
		filaments[0].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[0].length();
		stateArrayLengthFilled += 6*filaments[0].length();

		x_coord = LfcmBox_x/2 + 0.7*LL; // 0.1*10;
		y_coord = LfcmBox_y/2 + 0.5*LL;
		z_coord = LfcmBox_z/2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {0.,-1.,0. }; // Pointing towards tail
		HeadNormal = {0.,0.,-1. };
		filaments[1].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[1].length();
		stateArrayLengthFilled += 6*filaments[1].length();

	#elif TestInitialisation3

		float x_coord, y_coord, z_coord, LL;
		LL = 15*2.2;
		x_coord = LfcmBox_x/2 + 0.45*LL; // 0.1*10;
		y_coord = LfcmBox_y/2;
		z_coord = LfcmBox_z/2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {-1.,0.,0. }; // Pointing towards tail
		HeadNormal = {0.,-1.,0. };
		filaments[0].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[0].length();
		stateArrayLengthFilled += 6*filaments[0].length();

		x_coord = LfcmBox_x/2 + 0.55*LL; // 0.1*10;
		y_coord = LfcmBox_y/2 + 0.5*LL;
		z_coord = LfcmBox_z/2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {0.,-1.,0. }; // Pointing towards tail
		HeadNormal = {0.,0.,-1. };
		filaments[1].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[1].length();
		stateArrayLengthFilled += 6*filaments[1].length();

	#elif TestInitialisation4

		float x_coord, y_coord, z_coord, LL;
		LL = 100*2.2;
		x_coord = LfcmBox_x/2 - 0.5*LL; // 0.1*10;
		y_coord = LfcmBox_y/2;
		z_coord = LfcmBox_z/2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {1.,0.,0. }; // Pointing towards tail
		HeadNormal = {0.,1.,0. };
		filaments[0].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[0].length();
		stateArrayLengthFilled += 6*filaments[0].length();

		x_coord = LfcmBox_x/2; // 0.1*10;
		y_coord = LfcmBox_y/2 + 0.5*LL;
		z_coord = LfcmBox_z/2 + 2.2;
		HeadPosition = {x_coord, y_coord, z_coord};
		HeadTangent = {0.,-1.,0. }; // Pointing towards tail
		HeadNormal = {0.,0.,-1. };
		filaments[1].initialFilamentSetup(
			HeadPosition, HeadTangent, HeadNormal,
			forceArrayLengthFilled,  stateArrayLengthFilled);
		forceArrayLengthFilled += 1*filaments[1].length();
		stateArrayLengthFilled += 6*filaments[1].length();


	#else
		// Initial positions for swimmers.
		// continues placing swimmer until no overlap detected.
		for(int nn = 0; nn<Nsw; ++nn) {
			bool Overlap(true);
			do {
				HeadPosition = (InitialSeparation*nn)*HeadNormal;
				cout << "Proposing head position ... \n" << endl
				     << HeadPosition << endl << endl;
				filaments[nn].initialFilamentSetup(
					HeadPosition, HeadTangent, HeadNormal,
					forceArrayLengthFilled, stateArrayLengthFilled);
				Overlap = filaments[nn].checkOverlapWith(filaments,nn-1);
			}
			while(Overlap);
			cout << "... accepted head position."<<endl<<endl;
			forceArrayLengthFilled += 1*filaments[nn].length();
			stateArrayLengthFilled += 6*filaments[nn].length();
		}

		// ---------------------------------------------------------------------
	#endif

	return make_tuple(forceArrayLengthFilled, stateArrayLengthFilled);
}

// =============================================================================

tuple<int, int> two_bead_initialisation(vector<Filament>& filaments,
                                        double sep){

	// so that filament knows where it sits in the force vectors etc.
	int forceArrayLengthFilled = 0;
	// so that filament knows where it sits in the state vectors etc.
	int stateArrayLengthFilled = 0;
	float x_coord; float y_coord; float z_coord;
	int nn = 0;
	float kap_in; float c_in;

	vec HeadPosition(3);
	vec Head2Position(3);
	vec randVec(3);
	vec sepVec(3);

// Put two beads in first
	//arma_rng::set_seed(randSeed);
	bool TwoBeadsFoundSpace(false);
	do {
		// cout << "[Point T1]" << endl;
		// cout << "Length of filaments is " << filaments.size() << endl;
		// cout << "N_sw is " << N_sw << endl;
		Filament new_filament(1);
		filaments.push_back(new_filament);
		// cout << "Length of filaments is now " << filaments.size() << endl;
		// cout << "---" << endl;
		bool Overlap(true);

		cout << "Proposing CoM position of bead "
		     << "1/2"
		     << "..." << endl;
		do {
			// Uniform random in box of size RandomMessBoxSize
			HeadPosition.randu(); // Actually CoM position, not head pos'n
			HeadPosition[0] *= RandomMessBoxSize_x;
			HeadPosition[1] *= RandomMessBoxSize_y;
			HeadPosition[2] *= RandomMessBoxSize_z;

			// HeadPosition[0] = 3.3638346564e+00;
			// HeadPosition[1] = 1.5721890619e+01;
			// HeadPosition[2] = 5.7902473733e+01;

			filaments[N_sw-2].initialFilamentSetupIsoRand(
				HeadPosition,
				forceArrayLengthFilled, stateArrayLengthFilled, 0);
			Overlap = filaments[N_sw-2].checkOverlapWithPeriodic(filaments,
			                                                     N_sw-2, 1.05);
			// cout << "Overlap 1: " << Overlap << endl;

			//                                                                                                   was 2
		}
		while(Overlap);
		// TwoBeadsFoundSpace = true;

		Filament new_filament2(1);
		filaments.push_back(new_filament2);
		Overlap = true;
		cout << "Proposing CoM position of bead "
		     << "2/2 a distance " << sep << " away"
		     << "..." << endl;
		// bead 2
		int bead_2_attempts = 0;
		do {
			// Uniform random in box of size RandomMessBoxSize
			randVec.randu(); // just to get some random numbers
			double phi = randVec[0] * 2 * 3.1415926;
			double theta = randVec[1] * 3.1415926;
			sepVec[0] = sep*cos(phi)*sin(theta);
			sepVec[1] = sep*sin(phi)*sin(theta);
			sepVec[2] = sep*cos(theta);
			Head2Position = HeadPosition + sepVec;

			// Head2Position[0] = -1.4236706087e+00 + LfcmBox_x;
			// Head2Position[1] = 1.6853198149e+01;
			// Head2Position[2] = 6.1838199411e+01;

			filaments[N_sw-1].initialFilamentSetupIsoRand(
				Head2Position,
				forceArrayLengthFilled, stateArrayLengthFilled, 0);
			Overlap = filaments[N_sw-1].checkOverlapWithPeriodic(filaments,
			                                                     N_sw-2, 1.05);
			// cout << "Overlap 2: " << Overlap << endl;
			bead_2_attempts++;
			TwoBeadsFoundSpace = true;
			if(bead_2_attempts > 200) {
				cout << "Couldn't find space for bead 2" << endl;
				Overlap = false;
				TwoBeadsFoundSpace = false;
				// Remove last two filaments
				filaments.pop_back();
				filaments.pop_back();
			}
			//                                                                                                   was 2
		}
		while(Overlap);
	}
	while (!TwoBeadsFoundSpace);

	cout << " done." << endl;
	forceArrayLengthFilled += 1*filaments[N_sw-1].length();
	forceArrayLengthFilled += 1*filaments[N_sw-2].length();
	stateArrayLengthFilled += 6*filaments[N_sw-1].length();
	stateArrayLengthFilled += 6*filaments[N_sw-2].length();

	cout << "Placed at " << HeadPosition << " and " << Head2Position << endl;

	cout << " 2 beads initialised randomly." << endl;
	return make_tuple(forceArrayLengthFilled, stateArrayLengthFilled);
}


// =============================================================================

/**
    apply_filament_properties   Applies properties of filaments, as
                                given by the switches in config.hpp.
                                NOTE: This function alters the vector
                                `filaments` .
 */
void apply_filament_properties(vector<Filament>& filaments) {
	// Only run by rank 0

	#if GTinitialisation
		// Initialisation in a cylinder, see Gustavsson and Tornberg paper.
		double swim_phase = 0;
		for(int nn = 0; nn<Nsw; ++nn) {
			filaments[nn].Phase = swim_phase;
			swim_phase += PI/8+PI/2;
		}

		// ---------------------------------------------------------------------

	#elif RANDinitialisation
		// Random cloud of particles
		int randSeed = 2;
		arma_rng::set_seed(randSeed);

		double swim_phase = 0;
		for(int nn = 0; nn<Nsw; ++nn) {
			filaments[nn].Phase = swim_phase;
			vec tmprand(1,fill::zeros);
			tmprand.randu();
			swim_phase = 2*PI*tmprand[0];
		}
		// ---------------------------------------------------------------------

	#elif RandomMessToSwimThroughInitialisation || \
	RandomMessToSwimThroughInitialisationAllowOverlap
		// Random mess to swim
		double frozen_time = 0;

		// Put swimmer in first
		for(int nn = 0; nn<RandomMessNumSwimmingFilaments; ++nn) {
			filaments[nn].Phase = -PI/2; // Make them all have sin waves not cos waves
			filaments[nn].Alpha = SwimmingHelixAlpha;
			filaments[nn].Beta = SwimmingHelixBeta;
			filaments[nn].myKAP = Kap;
			filaments[nn].myC = C;

			if(FrozenSwimmer) {
				filaments[nn].Phase = -PI/8;
				filaments[nn].Phase -= filaments[nn].AngFrequency*frozen_time;
				filaments[nn].AngFrequency = 0;
			}
		}

		// Network
		for(int nn = RandomMessNumSwimmingFilaments; nn<N_sw; ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = 0;
			filaments[nn].Beta = 0;
			filaments[nn].myKAP = Kap_NetworkFilaments;
			filaments[nn].myC = C_NetworkFilaments;
		}
		// ---------------------------------------------------------------------


	#elif TwoParticleInitialisation
		// Two particle tests

		// Network
		for(int nn = 0; nn<N_sw; ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = 0;
			filaments[nn].Beta = 0;
			filaments[nn].myKAP = Kap_NetworkFilaments;
			filaments[nn].myC = C_NetworkFilaments;
		}
		// ---------------------------------------------------------------------

	#elif ConnectedNodesInitialisation || \
	ConnectedNodesInitialisationAllowOverlap
		// Random mess to swim

		// Put swimmer in first
		for(int nn = 0; nn<ConnectedNodesNumSwimmingFilaments; ++nn) {
			filaments[nn].Phase = -PI/2; // Make them all have sin waves not cos waves
			filaments[nn].Alpha = SwimmingHelixAlpha;
			filaments[nn].Beta = SwimmingHelixBeta;
			filaments[nn].myKAP = Kap;
			filaments[nn].myC = C;
		}

		// Network
		for(int nn = ConnectedNodesNumSwimmingFilaments; nn<filaments.size(); ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = 0;
			filaments[nn].Beta = 0;
			filaments[nn].myKAP = Kap_NetworkFilaments;
			filaments[nn].myC = C_NetworkFilaments;
		}
		// ---------------------------------------------------------------------

	#elif GTTallBoxInitialisation
		// Tall box for sedimenting filaments in periodic domain.
		// Inspired by G&T fig 5.

		// No swimming allowed
		for(int nn = 0; nn<N_sw; ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = 0;
			filaments[nn].Beta = 0;
			filaments[nn].myKAP = Kap;
			filaments[nn].myC = C;
		}

		// ---------------------------------------------------------------------
	#elif TempInitialisation
		// Random mess to swim

		// Put swimmer in first
		for(int nn = 0; nn<RandomMessNumSwimmingFilaments; ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = SwimmingHelixAlpha;
			filaments[nn].Beta = SwimmingHelixBeta;
			filaments[nn].myKAP = Kap;
			filaments[nn].myC = C;
		}

		// Network
		for(int nn = RandomMessNumSwimmingFilaments; nn<N_sw; ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = 0;
			filaments[nn].Beta = 0;
			filaments[nn].myKAP = Kap_NetworkFilaments;
			filaments[nn].myC = C_NetworkFilaments;
		}
		// ---------------------------------------------------------------------

	#elif StackedInitialisation
		// Sink
		for(int nn = 0; nn<N_sw; ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = 0;
			filaments[nn].Beta = 0;
			filaments[nn].myKAP = Kap;
			filaments[nn].myC = C;
		}

		// ---------------------------------------------------------------------

	#elif TestInitialisation1 || TestInitialisation2 || TestInitialisation3 || TestInitialisation4
		// Random mess to swim

		// Put swimmer in first
		for(int nn = 0; nn<RandomMessNumSwimmingFilaments; ++nn) {
			filaments[nn].Phase = -PI/2; // Make them all have sin waves not cos waves
			filaments[nn].Alpha = SwimmingHelixAlpha;
			filaments[nn].Beta = SwimmingHelixBeta;
			filaments[nn].myKAP = Kap;
			filaments[nn].myC = C;
		}

		// Network
		for(int nn = RandomMessNumSwimmingFilaments; nn<N_sw; ++nn) {
			filaments[nn].Phase = 0;
			filaments[nn].Alpha = 0;
			filaments[nn].Beta = 0;
			filaments[nn].myKAP = Kap_NetworkFilaments;
			filaments[nn].myC = C_NetworkFilaments;
		}
		// ---------------------------------------------------------------------

	#else
		// Random swim phase
		for(int nn = 0; nn<Nsw; ++nn) {
			vec tmprand(1,fill::zeros);
			tmprand.randu();
			filaments[nn].Phase = 2*PI*tmprand[0];
		}

	#endif
	cout << "Filament properties applied." << endl;
}

// =============================================================================
