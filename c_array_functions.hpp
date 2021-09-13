#ifndef C_ARRAY_INCLUDED
#define C_ARRAY_INCLUDED

#include "multi_filament_header.hpp"

double *double_array_1D(int nx);

double **double_array_2D(int nx, int ny);

double ***double_array_3D(int nx, int ny, int nz, int total_local_size);

int *int_array_1D(int nx);

int **int_array_2D(int nx, int ny);

int ***int_array_3D(int nx, int ny, int nz);

void free_double_array_1D(double *array);

void free_double_array_2D(double **array);

void free_double_array_3D(double ***array, int nx);

void free_int_array_1D(int *array);

void free_int_array_2D(int **array);

void free_int_array_3D(int ***array, int nx);

void make_zero_1D(double *A, int N);

void make_zero(double **A, int N, int M);

void make_zero_3D(double ***A, int N, int M, int Q);

// Linkedlist

void box(double **Y);

#endif
