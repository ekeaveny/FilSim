#ifndef C_ARRAY_INCLUDED
#define C_ARRAY_INCLUDED

#include "multi_filament_header.hpp"

double *double_array_1D(int nx){
	double *array;
	array = (double *) malloc(nx * sizeof(double));
	return array;
}

double **double_array_2D(int nx, int ny){
	int i = 0;
	double **array;
	array = (double **) malloc(nx * sizeof(double *));
	array[0] = (double *) malloc(nx*ny * sizeof(double));
	for(i = 1; i < nx; i++) {
		array[i] = array[0] + i * ny;
	}
	return array;
}

double ***double_array_3D(int nx, int ny, int nz, int total_local_size)
{
	int i,j;
	double ***arr;

	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// cout << "[Point $3a-XX1] " << myrank << ": " << nx*ny*nz << " " << total_local_size << endl;
	// allocate pointers to slices
	arr=(double ***) malloc(sizeof(double**)*nx);
	if (!arr) {
		cout << "Allocation error in double_array_3D(), code point 1" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	//cout << "[Point $3a-XX2] " << myrank << ": " << nx*ny*nz << " " << total_local_size << endl;
	// allocate rows and set pointers to them
	arr[0]=(double **) malloc(sizeof(double*)*nx*ny);
	if (!arr[0]) {
		cout << "Allocation error in double_array_3D(), code point 2" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	//cout << "[Point $3a-XX3] " << myrank << ": " << nx*ny*nz << " " << total_local_size << endl;
	//arr[0][0]=(double *) malloc(total_local_size*sizeof(double));
	arr[0][0]=(double *) calloc(total_local_size, sizeof(double));
	//arr[0][0]=(double *) calloc(nx*ny*nz, sizeof(double));
	if (!arr[0][0]) {
		cout << "Allocation error in double_array_3D(), code point 3" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	//cout << "[Point $3a-XX4] " << myrank << ": " << nx*ny*nz << " " << total_local_size << endl;
	// Note that these arrays won't allocate a contiguous region of memory.
	// So we add, so long as nx > 0:
	if (nx>0) {
		for(j=1; j<ny; j++) {
			arr[0][j] = arr[0][j-1] + nz;
		}
	}
	//cout << "[Point $3a-XX5] " << myrank << ": " << nx*ny*nz << " " << total_local_size << endl;
	for(i=1; i<nx; i++) {
		arr[i] = arr[i-1] + ny;
		arr[i][0] = arr[i-1][0] + ny*nz;
		for(j=0; j<ny; j++) {
			arr[i][j] = arr[i][j-1] + nz;
		}
	}
	//cout << "[Point $3a-XX6] " << myrank << ": " << nx*ny*nz << " " << total_local_size << endl;

	// return pointer
	return arr;
}

/*
   double ***double_array_3D(int nx, int ny, int nz){
    double ***array;
    array = (double ***) malloc(nx * sizeof(double **));
    for(int i = 0; i < nx; i++) {
        array[i] = (double **) malloc(ny * sizeof(double *));
    }
    if (nx*ny*nz > 0){ // In case we are initialising a zero array.
        array[0][0] = (double *) malloc(nx*ny*nz * sizeof(double));
    }
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            array[i][j] = array[0][0] + i*ny*nz + j*nz;
        }
    }
    return array;
   }*/

int *int_array_1D(int nx){
	int *array;
	array = (int *) malloc(nx * sizeof(int));
	return array;
}

int **int_array_2D(int nx, int ny){
	int i = 0;
	int **array;
	array = (int **) malloc(nx * sizeof(int *));
	array[0] = (int *) malloc(nx*ny * sizeof(int));
	for(i = 1; i < nx; i++) {
		array[i] = array[0] + i * ny;
	}
	return array;
}

int ***int_array_3D(int nx, int ny, int nz){
	int i = 0, j = 0;
	int ***array;
	array = (int ***) malloc(nx * sizeof(int **));
	for(i = 0; i < nx; i++) {
		array[i] = (int **) malloc(ny * sizeof(int *));
	}
	array[0][0] = (int *) malloc(nx*ny*nz * sizeof(int));
	for(i = 0; i < nx; i++) {
		for(j = 0; j < ny; j++) {
			array[i][j] = array[0][0] + i*ny*nz + j*nz;
		}
	}
	return array;
}

void free_double_array_1D(double *array){
	free(array);
	return;
}

void free_double_array_2D(double **array){
	free(array[0]);
	free(array);
	return;
}

/*
   void free_double_array_3D(double ***array, int nx){
    int i=0;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    cout << "freeing A " << myrank << " " << nx << endl;
    free(array[0][0]);
    cout << "freeing B " << myrank << " " << nx << endl;
    for(i = 0; i < nx; i++) {
        cout << "freeing Ba... " << myrank << " " << nx << " " << i << endl;
        free(array[i]);
        cout << "...done" << endl;
    }
    cout << "freeing C " << myrank << " " << nx << endl;
    free(array);
    cout << "freeing D " << myrank << " " << nx << endl;
    return;
   } */

// free a double array 3d as initialised in the new way.
// (from https://github.com/golosio/xrmc/blob/master/src/arrayNd/arrayNd.cpp)
void free_double_array_3D(double ***array, int nx){
	free(array[0][0]);
	free(array[0]);
	free(array);
	return;
}

void free_int_array_1D(int *array){
	free(array);
	return;
}

void free_int_array_2D(int **array){
	free(array[0]);
	free(array);
	return;
}

void free_int_array_3D(int ***array, int nx){
	int i=0;
	free(array[0][0]);
	for(i = 0; i < nx; i++) {
		free(array[i]);
	}
	free(array);
	return;
}

void make_zero_1D(double *A, int N){
	for(int i = 0; i < N; i++) {
		A[i] = 0.0;
	}
	return;
}

void make_zero(double **A, int N, int M){
	int i,j;
	for(i = 0; i < N; i++) {
		for(j = 0; j < M; j++) {
			A[i][j] = 0.0;
		}
	}
	return;
}

void make_zero_3D(double ***A, int N, int M, int Q){
	int i,j,k;
	for(i = 0; i < N; i++) {
		for(j = 0; j < M; j++) {
			for(k = 0; k < Q; k++) {
				A[i][j][k] = 0.0;
			}
		}
	}
	return;
}




// Linkedlist


/*
   int longest_filament(vector<Filament>& filaments){
    // Find longest filament
    max_Nw = 0;
    for(int i_fil=0; i_fil < Nsw; ++i_fil) {
        if (filaments[i_fil].myNworm > max_Nw) {
            max_Nw = filaments[i_fil].myNworm;
        }
    }
    return max_Nw;
   }
 */


void box(double **Y){
	double xr=0.0, yr=0.0, zr=0.0;
	int j = 0;
	double LX = 2.0*PI;
	for(j = 0; j<N; j++) {
		Y[j][0] = Y[j][0] - PI2 * floor( Y[j][0] /  PI2);
		Y[j][1] = Y[j][1] - PI2 * floor( Y[j][1] /  PI2);
		Y[j][2] = Y[j][2] - PI2 * floor( Y[j][2] /  PI2);
	}
	return;
}


#endif
