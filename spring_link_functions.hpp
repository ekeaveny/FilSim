// inclusion guards
#ifndef SPRING_LINK_FUNCTIONS_IS_INCLUDED
#define SPRING_LINK_FUNCTIONS_IS_INCLUDED

#include <math.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <armadillo>
//#include <vector>

#include "multi_filament_header.hpp"
//#include "UnitQuaternion.hpp"
//#include "Bead.hpp"
//#include "Filament.hpp"

using namespace std;
using namespace arma;

struct spring_links {

	vector< vector<int> > spring_links;

	vector<vector<vector<int> > > filaments_at_nodes;

	void read_spring_links_from_file(string data_file_name);

	void save_spring_links_to_file(string data_file_name);

	int size();

	void decide_spring_links(vector<Filament>& filaments);

	void apply_spring_link_forces(vector<Filament>& filaments);

	void set_num_nodes(int num_nodes);
};

// vector< vector<int> > decide_spring_links(vector<Filament>& filaments,
//                                           vector<vector<vector<int> > >& filaments_at_nodes);

// void apply_spring_link_forces(vector<Filament>& filaments,
//                               vector< vector<int> >& spring_links);

bool checkNodeOverlapWithPeriodic(vector<vec>& Nodes,
                                  const int NumberOthers,
                                  float separation_allowed);

vec VecPeriodicX(const vec& XX, const vec& YY);

vector< vector<int> > decide_node_to_node_filaments(vector<vec>& Nodes);

tuple<int,vec,vec> beads_between_nodes(vec node1, vec node2);

#endif
