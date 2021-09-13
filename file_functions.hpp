// inclusion guards
#ifndef FILE_FUNCTIONS_IS_INCLUDED
#define FILE_FUNCTIONS_IS_INCLUDED

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

void save_parameter_values_to_file(string parameter_file_name,
                                   double LfcmBox_x, double LfcmBox_y,
                                   double LfcmBox_z,
                                   ivec BeadNumbers, vec BendingFactors,
                                   vector<Filament>& filaments);

void save_data_column_names_to_file(string data_file_name,
                                    vector<Filament>& filaments);

void save_data_column_names_to_swimmer_velocity_file(string data_file_name,
                                                     vector<Filament>& filaments);

void save_data_column_names_to_swimmer_collision_file(string data_file_name,
                                                      vector<Filament>& filaments);

void save_data_column_names_to_forces_file(string data_file_name,
                                           vector<Filament>& filaments);

void save_data_to_file(string data_file_name,
                       vector<Filament>& filaments,
                       int nt);

void save_force_data_to_file(string data_file_name,
                             vector<Filament>& filaments,
                             double **F, double **T,
                             int nt);

void save_swimmer_velocity_data_to_file(string data_file_name,
                                        vector<Filament>& filaments,
                                        double **V,double **W,double **F,
                                        int nt);

void save_swimmer_collision_data_to_file(string data_file_name,
                                         vector<Filament>& filaments,
                                         double F[][3],
                                         int nt);

void save_backup_data_to_file(string data_file_name,
                              vector<Filament>& filaments,
                              int nt);

tuple<int, int, int> continue_from_checkpoint(string backup_data_file_name,
                                              vector<Filament>& filaments);

#endif
