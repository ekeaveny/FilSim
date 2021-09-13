// inclusion guards
#ifndef FILAMENT_INITIALISATION_FUNCTIONS_IS_INCLUDED
#define FILAMENT_INITIALISATION_FUNCTIONS_IS_INCLUDED

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


#if EnableSpringLinks
	tuple<int, int> filament_position_initialisation(vector<Filament>& filaments,
	                                                 spring_links& spring_links);
#else
	tuple<int, int> filament_position_initialisation(vector<Filament>& filaments);
#endif

tuple<int, int> two_bead_initialisation(vector<Filament>& filaments,
                                        double sep);

void apply_filament_properties(vector<Filament>& filaments);

#endif
