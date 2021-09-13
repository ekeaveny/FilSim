#include "spring_link_functions.hpp"
#include "config.hpp" // Include parameters set by config.hpp

// Decide which links should exist
void spring_links::decide_spring_links(vector<Filament>& filaments){
	vec disp;
	double dist;
	// Only called from rank 0.

	// Using filaments_at_nodes, create spring links between them.
	for (int node = 0; node < this->filaments_at_nodes.size(); ++node) {
		for ( int end1=0; end1 < this->filaments_at_nodes[node].size(); ++end1) {
			for (int end2=end1+1; end2 < this->filaments_at_nodes[node].size(); ++end2) {
				this->spring_links.push_back({this->filaments_at_nodes[node][end1][0],
				                              this->filaments_at_nodes[node][end1][1],
				                              this->filaments_at_nodes[node][end2][0],
				                              this->filaments_at_nodes[node][end2][1]});
			}
		}
	}
}

// Apply the spring forces to the links which have already been decided.
void spring_links::apply_spring_link_forces(vector<Filament>& filaments){

	// This can totally be optimised, for example by storing the separation
	// distance in spring_links.
	int nn, bead_n, mm, bead_m;
	double dist, spring_mag;
	vec disp, spring_force;
	for (int link_i = 0; link_i < spring_links.size(); ++link_i) {
		nn = this->spring_links[link_i][0];
		bead_n = this->spring_links[link_i][1];
		mm = this->spring_links[link_i][2];
		bead_m = this->spring_links[link_i][3];
		disp = VecPeriodicX(filaments[nn].beads[bead_n].Xs,
		                    filaments[mm].beads[bead_m].Xs);
		dist = norm(disp);
		spring_mag = NetworkSpringConstant*(dist-NetworkSpringNaturalLength);
		spring_force = -spring_mag*disp/norm(disp);
		filaments[nn].beads[bead_n].F += spring_force;
		filaments[mm].beads[bead_m].F -= spring_force;
		//cout << "Connecting " << nn << "," << bead_n << " with " << mm << "," << bead_m << " (dist " << dist << ")" << endl;
	}
}

// =============================================================================

int spring_links::size(){
	return this->spring_links.size();
}

// Write spring link data to file
void spring_links::save_spring_links_to_file(string data_file_name){
	// Only called by rank 0.

	ofstream OutputFile (data_file_name); // Overwrite, don't append.
	OutputFile << "Spring link data file. Number of links written first, "
	           << "then springs between {filament A, bead B} and "
	           << "{filament C, bead D} in the form 'A B C D'." << endl;
	OutputFile << this->size() << endl;
	for(int nn = 0; nn<this->size(); ++nn) {
		OutputFile << this->spring_links[nn][0] << " "
		           << this->spring_links[nn][1] << " "
		           << this->spring_links[nn][2] << " "
		           << this->spring_links[nn][3] << "\n";
	}
	OutputFile.close();
}

// Read spring link data from file
void spring_links::read_spring_links_from_file(string data_file_name){

	string spring_link_data_file_name =  OutputFolder + data_file_name + "-springlinks.dat";
	ifstream input_file (spring_link_data_file_name);
	try {
		string s, num_links_as_string, num_fils_as_string;
		int num_links;
		getline(input_file,s); // First \n-delimited item. Ignore
		getline(input_file,num_links_as_string); // 2nd line
		num_links = stoi(num_links_as_string);
		for (int nn = 0; nn < num_links; ++nn) {
			vector<int> spring_link;
			// Items 1 to 3 on the row:
			for (int mm = 0; mm < 3; ++mm) {
				getline(input_file,s,' ');
				spring_link.push_back(stoi(s));
			}
			// And then for the 4th: item on the row
			getline(input_file,s);
			spring_link.push_back(stoi(s));
			this->spring_links.push_back(spring_link);
		}
	}
	catch (const std::exception& e) {
		cout << endl << "ERROR: Spring link file '" << spring_link_data_file_name
		     << "' either does not exist or is corrupted." << endl << endl;
		exit(-1);
	}
	input_file.close();
}

void spring_links::set_num_nodes(int num_nodes){
	this->filaments_at_nodes.resize(num_nodes);
}

// =============================================================================

/**
    checkNodeOverlapWithPeriodic   Check all nodes in Nodes and see that they
                                   are far enough apart from each other.
 */
bool checkNodeOverlapWithPeriodic(vector<vec>& Nodes,
                                  const int NumberOthers,
                                  float separation_allowed)  {
	// Note `separation_allowed` is a distance; not a separation factor, as in
	// the `checkOverlapWithPeriodic` function within the Filament class.

	vec com1(3);
	vec com2(3);
	float distance_apart;

	for(int n = 0; n<NumberOthers; ++n) {
		com1 = Nodes[NumberOthers];
		com2 = Nodes[n];
		distance_apart = norm(VecPeriodicX(com1, com2));
		if(distance_apart < separation_allowed) {
			return true;
		}
	}
	return false;
}

// =============================================================================

/**
    VecPeriodicX   Periodic vector in direction X-Y;
 */
vec VecPeriodicX(const vec& XX, const vec& YY){ //
	vec XY = XX - YY;
	XY[0] = XY[0] - LfcmBox_x * floor(XY[0]/LfcmBox_x + 0.5);
	XY[1] = XY[1] - LfcmBox_y * floor(XY[1]/LfcmBox_y + 0.5);
	XY[2] = XY[2] - LfcmBox_z * floor(XY[2]/LfcmBox_z + 0.5);
	return XY;
}

// =============================================================================

/**
    decide_node_to_node_filaments     Decide which pairs of nodes should
                                      have filaments between them.
 */
vector< vector<int> > decide_node_to_node_filaments(vector<vec>& Nodes){
	vector< vector<int> > node_to_node_filaments;
	vec disp;
	for (int nn = 0; nn < Nodes.size(); ++nn) {
		for (int mm = nn+1; mm < Nodes.size(); ++mm) {
			double rand = randu<double>();
			if (rand < ConnectedNodesConnectionProbability){
				disp = VecPeriodicX(Nodes[nn], Nodes[mm]);
				if (norm(disp) < ConnectedNodesConnectionRadius){			
					node_to_node_filaments.push_back({nn,mm});
				}
			}
		}
	}
	return node_to_node_filaments;
}

tuple<int,vec,vec> beads_between_nodes(vec node1, vec node2){
	vec node_to_node = VecPeriodicX(node2,node1); // node2-node1
	double node_to_node_dist = norm(node_to_node);
	int num_beads_required = (node_to_node_dist -
	                          2*ConnectedNodesSpacingAwayFromNode) / dL; // (Note this rounds down)
	vec CoM = node1 + node_to_node/2;
	vec unit_vector = node_to_node/node_to_node_dist; // node2-node1
	return make_tuple(num_beads_required, CoM, unit_vector);
}
