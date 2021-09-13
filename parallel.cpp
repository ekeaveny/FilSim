#include <iostream>
#include <stdlib.h>
#include <mpi.h>
using namespace std;



// from Eric's C code, but works in C++
double ***double_array_3D(int nx, int ny, int nz){
	int i = 0, j = 0;
	double ***array;
	array = (double ***) malloc(nx * sizeof(double **));
	for(i = 0; i < nx; i++) {
		array[i] = (double **) malloc(ny * sizeof(double *));
	}
	array[0][0] = (double *) malloc(nx*ny*nz * sizeof(double));
	for(i = 0; i < nx; i++) {
		for(j = 0; j < ny; j++) {
			array[i][j] = array[0][0] + i*ny*nz + j*nz;
		}
	}
	return array;
}

// from Eric's C code, but works in C++
void free_double_array_3D(double ***array, int nx){
	int i=0;
	free(array[0][0]);
	for(i = 0; i < nx; i++) {
		free(array[i]);
	}
	free(array);
	return;
}




int main(int argc, char **argv)
{
	int totalnodes, myrank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// std::cout << argv << endl;

	for(int i = 1; i< argc; ++i) {
		std::cout << argv[i] << endl;
	}



	if(myrank==0) {
		std::cout << "my rank is 0"<< endl;
	}
	else {
		std::cout << "my non-zero rank is " << myrank << endl;
	}

	int N = 5;

	double ***Y;
	Y = double_array_3D(N, 3, 4);

	if(myrank==0) {
		for(int i = 0; i<N; ++i) {
			for(int j = 0; j<3; ++j) {
				for(int k = 0; k<3; ++k) {
					Y[i][j][k] = i*j*k;
				}
			}
		}

		double serial_sum = 0;
		for(int i = 0; i<N; ++i) {
			for(int j = 0; j<3; ++j) {
				for(int k = 0; k<3; ++k) {
					serial_sum += Y[i][j][k];
				}
			}
		}
		std::cout << "serially computed sum = " << serial_sum <<endl;

	}

	// int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root,
	//              MPI_Comm comm )
	MPI_Bcast(&Y[0][0][0], 3*4*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	std::cout << "element " << myrank*Y[1][1][1] << endl;

	double par_sum = 0;
	for(int i = 0; i<N; ++i) {
		for(int j = 0; j<3; ++j) {
			par_sum += Y[i][j][myrank];
		}
	}

	std::cout<< "my partial sum is " << par_sum << endl;

	double total_sum;

	// all to all
	MPI_Allreduce(&par_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// all to root/master
	MPI_Reduce(&par_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	if(myrank == 0) {
		std::cout<< "total sum is " << total_sum << endl;
	}


	free_double_array_3D(Y, N);

	MPI_Finalize();
	return 0;
}
