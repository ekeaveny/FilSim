// inclusion guards
#ifndef PRINT_FUNCTIONS_IS_INCLUDED
#define PRINT_FUNCTIONS_IS_INCLUDED

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

string format_time(float elapsed_time);

void print_parameter_values_in_table(int Np, double LfcmBox_x, double LfcmBox_y,
                                     double LfcmBox_z, double L,
                                     double omega, double dt);

void set_precision_of_screen_output(int precision);

#endif
