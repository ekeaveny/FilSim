#ifndef FCM_INCLUDED
#define FCM_INCLUDED

#include "multi_filament_header.hpp"
#include "profilers.hpp"

struct Fcm {
	// The main guys here are: Y (current positions)
	//						   F (current forces)
	//						   T (current torques)

	/* Particle information - current, past and temporary. */
	double **Y, **Yi, **YD, **YDi, **V, **F, **T, **F0, **T0;
	double **VTEMP, **FTEMP, **V0, **W, **W0, **WTEMP, **p;
	double **AVG_TAU, **AVG_TAUTEMP;
	double **GA, **GS0, **Gp, **ETEMP, **E0, **E;

	/* Domain variables - velocity field, force distribution*/
	double ***ux, ***uy, ***uz, ***fx, ***fy, ***fz;
	//double *fx2, *fy2, *fz2;

	/* FFT associated variables - plans, wave vectors, local process information, padding */
	double ***workspace;
	double *qx, *qy, *qpad, *qxsq, *qysq, *qpadsq;
	int local_nx, local_x_start, local_nyt, local_y_start, total_local_size, pad;
	long total_points;
	rfftwnd_mpi_plan plan, iplan;

	// /* Grid spacing */
	// double dx, x, y, z;

	/* Variables related to the gaussian distributions associated with FCM and the magnetic dipole calculations */
	int ngdh;
	double **gaussx, **gaussy, **gaussz;
	double **gaussx_dip, **gaussy_dip, **gaussz_dip;
	double **grad_gaussx_dip, **grad_gaussy_dip, **grad_gaussz_dip;

	int **indx, **indy, **indz;
	double fdipdim;

	void assign_filament_data(const std::vector<Filament>& filaments);

	void initialise(int myrank);

	#if PROFILING
		void mobility_solve(profiler *profiler);
	#else
		void mobility_solve();
	#endif

	void friction_mobility_solve();

	void free_memory();

	void print_data(string begin_string = "", string end_string = "", bool pause_after = false);

	void mpi_broadcast();

};

void GA_setup(double **GA, double **T);

void gaussian_setup(int ngdh, double **Y,
                    double **gaussx, double **gaussy, double **gaussz,
                    double **gaussx_dip, double **gaussy_dip, double **gaussz_dip,
                    double **grad_gaussx_dip, double **grad_gaussy_dip, double **grad_gaussz_dip, double sigmadipsq,
                    int **indx, int **indy, int **indz,
                    double anorm, double anorm2, double anormdip, double anormdip2, double dx);

void mono_and_dipole_force_distribution(double ***fx, double ***fy, double ***fz,
                                        double **GS, double **GA, double **F,
                                        double **gaussx, double **gaussy, double **gaussz,
                                        double **gaussx_dip, double **gaussy_dip, double **gaussz_dip,
                                        double **grad_gaussx_dip, double **grad_gaussy_dip, double **grad_gaussz_dip,
                                        int **indx, int **indy, int **indz,
                                        int local_x_start, int local_nx);

void flow_solve(double ***fx, double ***fy, double ***fz,
                double ***ux, double ***uy, double ***uz,
                double *qx, double *qy, double *qpad,
                double *qxsq, double *qysq, double *qpadsq,
                int local_nyt, int local_y_start, int pad);

void particle_velocities(double ***ux, double ***uy, double ***uz, double **VTEMP,
                         double **gaussx, double **gaussy, double **gaussz,
                         int **indx, int **indy, int **indz,
                         int local_x_start, int local_nx, double dx);

void particle_rotations(double ***ux, double ***uy, double ***uz, double **WTEMP,
                        double **gaussx_dip, double **gaussy_dip, double **gaussz_dip,
                        double **grad_gaussx_dip, double **grad_gaussy_dip, double **grad_gaussz_dip,
                        int **indx, int **indy, int **indz,
                        int local_x_start, int local_nx, double dx);

void FCMStokeslet(double **Y, double **F,double **T,double **V,double **W);

#endif
