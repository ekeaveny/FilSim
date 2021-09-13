#include "multi_filament_header.hpp"
#include "c_array_functions.hpp"
#include "profilers.hpp"

// Initialise variables which are used in the FCM solve. Most are behind-the-
// scenes but the main ones which are interfaced externally are
// fcm.F, fcm.T, (input forces and torques), and
// fcm.V and fcm.W, (output velocities and angular velocities) for the mobility
// solve.
void Fcm::initialise(int myrank){
	if (myrank == 0) {
		cout << "Setting up FCM/FFT grids... " << std::flush;
	}

	this->ngdh = ngd/2;
	this->plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, NPTS_X, NPTS_Y, NPTS_Z,
	                                     FFTW_REAL_TO_COMPLEX, FFTW_MEASURE);
	this->iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, NPTS_X, NPTS_Y, NPTS_Z,
	                                      FFTW_COMPLEX_TO_REAL, FFTW_MEASURE);
	rfftwnd_mpi_local_sizes(this->plan,
	                        &this->local_nx, &this->local_x_start,
	                        &this->local_nyt, &this->local_y_start,
	                        &this->total_local_size);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $1]" << endl;

	this->pad = 2 * (NPTS_Z/2 + 1);
	this->total_points = NPTS_X*NPTS_Y*NPTS_Z;

	// Memory allocation for pointer variables EXCEPT linked-lists
	this->Y = double_array_2D(Np, 3);
	this->Yi = double_array_2D(Np, 3);
	this->YD = double_array_2D(Np, 3);
	this->YDi = double_array_2D(Np, 3);
	this->F = double_array_2D(Np, 3);
	this->T = double_array_2D(Np, 3);
	this->F0 = double_array_2D(Np, 3);
	this->T0 = double_array_2D(Np, 3);
	this->V = double_array_2D(Np, 3);
	this->VTEMP = double_array_2D(Np, 3);
	this->FTEMP = double_array_2D(Np, 3);
	this->WTEMP = double_array_2D(Np, 3);
	this->V0 = double_array_2D(Np, 3);
	this->W0 = double_array_2D(Np, 3);
	this->W = double_array_2D(Np, 3);
	this->Gp = double_array_2D(Np, 6);
	this->GS0 = double_array_2D(Np, 6);
	this->GA = double_array_2D(Np, 6);
	this->E = double_array_2D(Np, 6);
	this->E0 = double_array_2D(Np, 6);
	this->ETEMP = double_array_2D(Np, 6);
	this->AVG_TAU = double_array_2D(Np, 6);
	this->AVG_TAUTEMP = double_array_2D(Np, 6);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << endl;
	//cout << "[Point $2]" << endl;

	// Create FFT plans for field & flow solvers and obtain local process
	// information for the FFT. It is assumed that NPTS is divisible by the
	// number of processors used (totalnodes).
	this->ux=double_array_3D(this->local_nx, NPTS_Y, this->pad, this->total_local_size);
	//cout << "[Point $2az]" << myrank << endl;
	//sleep(1);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	//cout << endl;
	//cout << "[Point $3a] " << myrank << ": " << this->local_nx << " " << NPTS_Y << " " << this->pad << " " << this->total_local_size <<  endl;
	// If you ask for 15x24 (i.e. where it doesn't divide NPTS^3), you end up
	// with an error that happens somewhere between HERE...
	this->uy=double_array_3D(this->local_nx, NPTS_Y, this->pad, this->total_local_size);
	//cout << "[Point $3az]" << myrank << endl;
	//sleep(1);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	//cout << "[Point $3b]" << endl;
	this->uz=double_array_3D(this->local_nx, NPTS_Y, this->pad, this->total_local_size);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $3c]" << endl;

	this->fx=double_array_3D(this->local_nx, NPTS_Y, this->pad, this->total_local_size);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $3d]" << endl;
	this->fy=double_array_3D(this->local_nx, NPTS_Y, this->pad, this->total_local_size);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $3e]" << endl;
	this->fz=double_array_3D(this->local_nx, NPTS_Y, this->pad, this->total_local_size);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	//this->fx2 = (double*) malloc(sizeof(double) * this->total_local_size);
	//this->fy2 = (double*) malloc(sizeof(double) * this->total_local_size);
	//this->fz2 = (double*) malloc(sizeof(double) * this->total_local_size);

	//cout << "[Point $4]" << endl;
	// ... and HERE!
	this->workspace=double_array_3D(this->local_nx, NPTS_Y, this->pad, this->total_local_size);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $5]" << endl;
	this->gaussx=double_array_2D(Np, ngd);
	this->gaussy=double_array_2D(Np, ngd);
	this->gaussz=double_array_2D(Np, ngd);
	this->gaussx_dip=double_array_2D(Np, ngd);
	this->gaussy_dip=double_array_2D(Np, ngd);
	this->gaussz_dip=double_array_2D(Np, ngd);
	this->grad_gaussx_dip=double_array_2D(Np, ngd);
	this->grad_gaussy_dip=double_array_2D(Np, ngd);
	this->grad_gaussz_dip=double_array_2D(Np, ngd);
	this->indx=int_array_2D(Np, ngd);
	this->indy=int_array_2D(Np, ngd);
	this->indz=int_array_2D(Np, ngd);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $6]" << endl;

	this->qx=double_array_1D(NPTS_X); // Wavenumber in x-direction
	this->qxsq=double_array_1D(NPTS_X); // ,, squared
	this->qy=double_array_1D(NPTS_Y); // Wavenumber in x-direction
	this->qysq=double_array_1D(NPTS_Y); // ,, squared

	this->qpad=double_array_1D(this->pad); // Wavenumber in z-direction
	this->qpadsq=double_array_1D(this->pad); // ,, squared

	this->fdipdim = 20.0 * PI * a*a*a / 3.0;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $7]" << endl;

	// Padding needed for the real/complex inplace transform is computed.
	double FourierFactor_x = TwoPI/LfcmBox_x; // this is 2*pi/L, in Matlab qxfac
	double FourierFactor_y = TwoPI/LfcmBox_y; // this is 2*pi/L, in Matlab qxfac
	double FourierFactor_z = TwoPI/LfcmBox_z; // this is 2*pi/L, in Matlab qxfac

	int npts_half =  (NPTS_X / 2);
	for(int i = 0; i < NPTS_X; i++) {
		if(i <= npts_half) {
			this->qx[i] = i*FourierFactor_x;
		} else {
			this->qx[i] = (i - NPTS_X)*FourierFactor_x;
		}
		this->qxsq[i] = this->qx[i]*this->qx[i];
	}
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $8]" << endl;

	npts_half =  (NPTS_Y / 2);
	for(int i = 0; i < NPTS_Y; i++) {
		if(i <= npts_half) {
			this->qy[i] = i*FourierFactor_y;
		} else {
			this->qy[i] = (i - NPTS_Y)*FourierFactor_y;
		}
		this->qysq[i] = this->qy[i]*this->qy[i];
	}
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point $9]" << endl;

	for(int i = 0; i < this->pad; i = i+2) {
		this->qpad[i] = ((double) (i/2))*FourierFactor_z;
		this->qpad[i+1] = this->qpad[i];
		this->qpadsq[i] = this->qpad[i] * this->qpad[i];
		this->qpadsq[i+1] = this->qpadsq[i];
	}

	MPI_Barrier(MPI_COMM_WORLD); // Only say 'done' if every thread is happy.
	if (myrank == 0) {
		cout << " done." << endl;
	}
};


// Solve mobility problem. For given fcm.F and fcm.T, find fcm.V_fcm and
// fcm.W_fcm (velocities and angular velocities). Requires initialising the fcm
// struct first.
#if PROFILING
void Fcm::mobility_solve(profiler *prof)
#else
void Fcm::mobility_solve()
#endif
{
	#if PROFILING
		//profiler profiler = *prof;
		int myrank;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		if(myrank == 0) {
			prof->fcm_make_zero.start();
			//auto elapsed_time = (std::chrono::high_resolution_clock::now() - profiler.totalTime.start_time);
			//int elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count();
			//cout << "Value " << elapsed_ms << endl;
		}
	#endif

	make_zero(this->V, Np, 3);
	make_zero(this->W, Np, 3);
	make_zero(this->VTEMP, Np, 3);
	make_zero(this->WTEMP, Np, 3);
	// make_zero(GA, N, 6);

	// dipole terms (only rotlets, no stresslets)
	GA_setup(this->GA, this->T);
	make_zero(this->GS0, Np, 6);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	//if (this->local_nx*NPTS_Y*this->pad > 0){

	make_zero_3D(this->fx, this->local_nx, NPTS_Y, this->pad);
	make_zero_3D(this->fy, this->local_nx, NPTS_Y, this->pad);
	make_zero_3D(this->fz, this->local_nx, NPTS_Y, this->pad);
	//make_zero_1D(this->fx2, this->total_local_size);
	//make_zero_1D(this->fy2, this->total_local_size);
	//make_zero_1D(this->fz2, this->total_local_size);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point Ka]" << myrank << endl;

	#if PROFILING
		if(myrank == 0) {
			prof->fcm_make_zero.end();
			prof->fcm_gaussian_setup.start();
		}
	#endif

	//  Set up gaussians for FCM
	gaussian_setup(this->ngdh, this->Y,
	               this->gaussx, this->gaussy, this->gaussz,
	               this->gaussx_dip, this->gaussy_dip, this->gaussz_dip,
	               this->grad_gaussx_dip, this->grad_gaussy_dip,
	               this->grad_gaussz_dip, sigmadipsq,
	               this->indx, this->indy, this->indz,
	               anorm, anorm2, anormdip, anormdip2, dx);

	#if PROFILING
		if(myrank == 0) {
			prof->fcm_gaussian_setup.end();
			//cout << prof->fcm_gaussian_setup.count << endl;
			prof->fcm_force_distribution.start();
		}
	#endif

	// Set up the FCM force distribution
	mono_and_dipole_force_distribution(
		this->fx, this->fy, this->fz, this->GS0, this->GA, this->F,
		this->gaussx, this->gaussy, this->gaussz,
		this->gaussx_dip, this->gaussy_dip, this->gaussz_dip,
		this->grad_gaussx_dip, this->grad_gaussy_dip, this->grad_gaussz_dip,
		this->indx, this->indy, this->indz, this->local_x_start, this->local_nx);//,
	//this->fx2, this->fy2, this->fz2);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point Kb]" << myrank << endl;

	/*
	   // TEST 10/7/19
	   double LLX = 2.0*PI/(double)NPTS_X;
	   double LLY = 2.0*PI/(double)NPTS_Y;
	   double LLZ = 2.0*PI/(double)NPTS_Z;
	   cout << LLX << " " << LLY << " " << LLZ << "!" << endl;
	   cin.get();
	   for (int ii=0; ii<this->local_nx; ++ii) {
	    for (int jj=0; jj<NPTS_Y; ++jj) {
	        for (int kk=0; kk<NPTS_Z; ++kk) {
	            //this->fx[ii][jj][kk] = (LLX*LLX + LLY*LLY + LLZ*LLZ)*std::sin(ii*LLX)*std::sin(jj*LLY)*std::sin(kk*LLZ);
	            this->fx[ii][jj][kk] = 0;
	            this->fy[ii][jj][kk] = 0; //std::sin(ii*LL)*std::sin(jj*LL)*std::sin(kk*LL);
	            this->fz[ii][jj][kk] = std::sin((double)(this->local_x_start+ii)*LLX); //std::sin(ii*LL)*std::sin(jj*LL)*std::sin(kk*LL);


	            //  cout << ii << "," << jj << "," << kk << ": "
	            //   << this->fx[ii][jj][kk] << ","
	            //   << this->fy[ii][jj][kk] << ","
	            //   << this->fz[ii][jj][kk] << endl;

	        }
	    }
	   }
	   /*
	   for (int ii=0; ii<this->local_nx; ++ii) {
	    for (int jj=0; jj<NPTS_XY; ++jj) {
	        for (int kk=0; kk<this->pad; ++kk) {
	            if (abs(fx[ii][jj][kk]) > 1e-10) {
	                cout << "f " << ii << " " << jj << " " << kk << ": " << fx[ii][jj][kk] << endl;
	            }
	        }
	    }
	   }
	 */

	//cout << "[Points Kbm]" << myrank << " " << &this->fx[0][0][0] << " // " << &this->fx[0][0] << " // " << &this->fx[0] << " // " << this->fx << endl;
	//cout << "[Points Kbn]" << myrank << " " << &this->ux[0][0][0] << " // " << &this->ux[0][0] << " // " << &this->ux[0] << " // " << this->ux << endl;
	//cout << "[Points Kbo]" << myrank << " " << &this->workspace[0][0][0] << " // " << &this->workspace[0][0] << " // " << &this->workspace[0] << " // " << this->workspace << endl;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	// FFT the FCM force distribution

	/*
	   double *p = &this->fz[0][0][0];
	   for (int x = 0; x < this->local_nx; ++x){
	    for (int y = 0; y < NPTS_Y; ++y){
	        for (int z = 0; z < NPTS_Z; ++z){
	            cout << "*p " << "[" << x << "][" << y << "][" << z << "] "<< *p++ << " " << fz[x][y][z] << endl;
	        }
	 * p++;
	 * p++;
	    }
	   }
	   cout << "Forty more ";
	   for (int zzz=0; zzz<40; ++zzz) cout << *p++ << " ";
	   cout << endl;
	   MPI_Barrier(MPI_COMM_WORLD);
	 */

	//double *data, *work;
	//data = (double*) malloc(sizeof(double) * this->total_local_size);
	//cout << "[Total local size]" << myrank << " " << this->total_local_size << endl;
	/* workspace is the same size as the data: */
	//work = (fftw_real*) malloc(sizeof(fftw_real) * this->total_local_size);
	/* initialize data to f(x,y,z): */
	/*
	   for (int x = 0; x < this->local_nx; ++x)
	        for (int y = 0; y < NPTS_Y; ++y)
	                for (int z = 0; z < NPTS_Z; ++z)
	                        data[(x*NPTS_Y + y) * this->pad + z]
	                                = fx[x][y][z];
	   // Now, compute the forward transform:
	   rfftwnd_mpi(this->plan, 1, data, &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	   //cout << "[Point Kby1]" << myrank << endl;
	   //MPI_Barrier(MPI_COMM_WORLD); // remove later
	   for (int x = 0; x < this->local_nx; ++x)
	        for (int y = 0; y < NPTS_Y; ++y)
	                for (int z = 0; z < NPTS_Z; ++z)
	                        data[(x*NPTS_Y + y) * this->pad + z]
	                                = fy[x][y][z];
	   // Now, compute the forward transform:
	   rfftwnd_mpi(this->plan, 1, data, &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	   cout << "[Point Kby2]" << myrank << endl;
	   MPI_Barrier(MPI_COMM_WORLD); // remove later
	   for (int x = 0; x < this->local_nx; ++x)
	        for (int y = 0; y < NPTS_Y; ++y)
	                for (int z = 0; z < NPTS_Z; ++z)
	                        data[(x*NPTS_Y + y) * this->pad + z]
	                                = fz[x][y][z];
	   // Now, compute the forward transform:
	   rfftwnd_mpi(this->plan, 1, data, &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	   cout << "[Point Kby3]" << myrank << endl;
	   MPI_Barrier(MPI_COMM_WORLD); // remove later
	 */





	#if PROFILING
		if(myrank == 0) {
			prof->fcm_force_distribution.end();
			prof->fcm_fft_forward.start();
		}
	#endif



	rfftwnd_mpi(this->plan, 1, &this->fx[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	//cout << "[Point Kbz1]" << myrank << endl;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	// What if I comment out the rest?

	//rfftwnd_mpi(this->plan, 1, &this->fy[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	rfftwnd_mpi(this->plan, 1, &this->fy[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	//cout << "[Point Kbz2]" << myrank << endl;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	//rfftwnd_mpi(this->plan, 1, &this->fz[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	rfftwnd_mpi(this->plan, 1, &this->fz[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	//cout << "[Point Kbz3]" << myrank << endl;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	#if PROFILING
		if(myrank == 0) {
			prof->fcm_fft_forward.end();
			prof->fcm_flow_solve.start();
		}
	#endif

	//cout << "XA" << endl;
	//cout << "[Point Kc]" << myrank << endl;
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	// Solve for the flow field
	flow_solve(this->fx, this->fy, this->fz,
	           this->ux, this->uy, this->uz,
	           this->qx, this->qy, this->qpad,
	           this->qxsq, this->qysq, this->qpadsq,
	           this->local_nyt, this->local_y_start, this->pad);
	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point Kca]" << myrank << endl;

	#if PROFILING
		if(myrank == 0) {
			prof->fcm_flow_solve.end();
			prof->fcm_fft_backward.start();
		}
	#endif

	// FFT the velocity field
	rfftwnd_mpi(this->iplan, 1, &this->ux[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	//cout << "[Point Kcb]" << myrank << endl;
	rfftwnd_mpi(this->iplan, 1, &this->uy[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);
	rfftwnd_mpi(this->iplan, 1, &this->uz[0][0][0], &this->workspace[0][0][0], FFTW_TRANSPOSED_ORDER);

	//MPI_Barrier(MPI_COMM_WORLD); // remove later
	//cout << "[Point Kd]" << myrank << endl;

	/*
	   // TEST
	   for (int ii=0; ii<this->local_nx; ++ii) {
	    for (int jj=0; jj<NPTS_Y; ++jj) {
	        for (int kk=0; kk<this->pad; ++kk) {
	            if ((abs(this->ux[ii][jj][kk]) > 1e-16 ||
	                 abs(this->uy[ii][jj][kk]) > 1e-16 ||
	                 abs(this->uz[ii][jj][kk]) > 1e-16)
	                && kk==0 // NPTS_X-1
	                && jj==0 //jj==NPTS_Y-1
	                ) {
	                cout << "u " << local_x_start+ii << "," << jj << "," << kk << ": "
	                     << this->ux[ii][jj][kk] << ","
	                     << this->uy[ii][jj][kk] << ","
	                     << this->uz[ii][jj][kk] << endl;
	            }
	        }
	    }
	   }
	 */

	#if PROFILING
		if(myrank == 0) {
			prof->fcm_fft_backward.end();
			prof->fcm_particle_velocities_rotations.start();
		}
	#endif

	// Compute the particle velocities.
	particle_velocities(this->ux, this->uy, this->uz, this->VTEMP,
	                    this->gaussx, this->gaussy, this->gaussz,
	                    this->indx, this->indy, this->indz,
	                    this->local_x_start, this->local_nx, dx);

	particle_rotations(this->ux, this->uy, this->uz, this->WTEMP,
	                   this->gaussx_dip, this->gaussy_dip, this->gaussz_dip,
	                   this->grad_gaussx_dip, this->grad_gaussy_dip,
	                   this->grad_gaussz_dip,
	                   this->indx, this->indy, this->indz,
	                   this->local_x_start, this->local_nx, dx);

	#if PROFILING
		if(myrank == 0) {
			prof->fcm_particle_velocities_rotations.end();
		}
	#endif
	//}
	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	MPI_Reduce(&this->VTEMP[0][0], &this->V[0][0], 3*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&this->WTEMP[0][0], &this->W[0][0], 3*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


};


// Free memory associated with all of the FCM variables.
void Fcm::free_memory(){
	rfftwnd_mpi_destroy_plan(this->plan);
	rfftwnd_mpi_destroy_plan(this->iplan);

	free_double_array_2D(this->Y);
	free_double_array_2D(this->Yi);
	free_double_array_2D(this->V);
	free_double_array_2D(this->F);
	free_double_array_2D(this->T);
	free_double_array_2D(this->VTEMP);
	free_double_array_2D(this->V0);
	free_double_array_2D(this->W);
	free_double_array_2D(this->WTEMP);
	free_double_array_2D(this->W0);
	free_double_array_2D(this->FTEMP);
	free_double_array_2D(this->E);
	free_double_array_2D(this->E0);
	free_double_array_2D(this->Gp);
	free_double_array_2D(this->ETEMP);
	free_double_array_2D(this->GS0);
	free_double_array_2D(this->GA);

	free_double_array_3D(this->ux, this->local_nx);
	free_double_array_3D(this->uy, this->local_nx);
	free_double_array_3D(this->uz, this->local_nx);
	free_double_array_3D(this->fx, this->local_nx);
	free_double_array_3D(this->fy, this->local_nx);
	free_double_array_3D(this->fz, this->local_nx);

	free_double_array_3D(this->workspace, this->local_nx);

	free_double_array_2D(this->gaussx);
	free_double_array_2D(this->gaussy);
	free_double_array_2D(this->gaussz);
	free_double_array_2D(this->gaussx_dip);
	free_double_array_2D(this->gaussy_dip);
	free_double_array_2D(this->gaussz_dip);
	free_double_array_2D(this->grad_gaussx_dip);
	free_double_array_2D(this->grad_gaussy_dip);
	free_double_array_2D(this->grad_gaussz_dip);

	free_int_array_2D(this->indx);
	free_int_array_2D(this->indy);
	free_int_array_2D(this->indz);

	free_double_array_1D(this->qx);
	free_double_array_1D(this->qy);
	free_double_array_1D(this->qpad);
	free_double_array_1D(this->qxsq);
	free_double_array_1D(this->qysq);
	free_double_array_1D(this->qpadsq);
}

void Fcm::print_data(string begin_string, string end_string, bool pause_after){
	cout << begin_string << endl;
	for(int np = 0; np < min(Np,200); ++np) {
		std::cout << "Y[" << np << "][0]: "<< this->Y[np][0]<<", Y[" << np << "][1]: "<< this->Y[np][1]<< ", Y[" << np << "][2]: "<< this->Y[np][2] << endl;
		std::cout << "F[" << np << "][0]: "<< this->F[np][0]<<", F[" << np << "][1]: "<< this->F[np][1]<< ", F[" << np << "][2]: "<< this->F[np][2] << endl;
		//std::cout << "T[" << np << "][0]: "<< this->T[np][0]<<", T[" << np << "][1]: "<< this->T[np][1]<< ", T[" << np << "][2]: "<< this->T[np][2] << endl;
		std::cout << "V[" << np << "][0]: "<< this->V[np][0]<<", V[" << np << "][1]: "<< this->V[np][1]<< ", V[" << np << "][2]: "<< this->V[np][2] << endl;
		//std::cout << "W[" << np << "][0]: "<< this->W[np][0]<<", W[" << np << "][1]: "<< this->W[np][1]<< ", W[" << np << "][2]: "<< this->W[np][2] << endl;
	}
	cout << end_string << endl;
	if (pause_after) {
		std::cin.get();
	}
}

void GA_setup(double **GA, double **T){
	int i;
	for(i = 0; i < Np; i++) {
		GA[i][0] = 0.0;
		GA[i][1] = 0.0;
		GA[i][2] = 0.0;
		GA[i][3] = 0.5*T[i][2];
		GA[i][4] = -0.5*T[i][1];
		GA[i][5] = 0.5*T[i][0];
	}
	return;
}

void gaussian_setup(
	int ngdh, double **Y,
	double **gaussx, double **gaussy, double **gaussz,
	double **gaussx_dip, double **gaussy_dip, double **gaussz_dip,
	double **grad_gaussx_dip, double **grad_gaussy_dip,
	double **grad_gaussz_dip, double sigmadipsq,
	int **indx, int **indy, int **indz,
	double anorm, double anorm2, double anormdip, double anormdip2,
	double dx){

	int np, i, xc, yc, zc;
	int xg, yg, zg;

	double xx;

	for(np = 0; np < Np; np++) {
		xc = (int) (Y[np][0]/dx);
		yc = (int) (Y[np][1]/dx);
		zc = (int) (Y[np][2]/dx);

		for(i = 0; i < ngd; i++) {
			xg = xc - ngdh + (i+1);
			indx[np][i] = xg - NPTS_X * ((int) std::floor( ((double) xg) / ((double) NPTS_X)));
			xx = ((double) xg)*dx-Y[np][0];
			gaussx[np][i] = anorm*std::exp(-xx*xx/anorm2);
			gaussx_dip[np][i] = anormdip*std::exp(-xx*xx/anormdip2);
			grad_gaussx_dip[np][i] = -xx / sigmadipsq;

			yg = yc - ngdh + (i+1);
			indy[np][i] = yg - NPTS_Y * ((int) std::floor( ((double) yg) / ((double) NPTS_Y)));
			xx = ((double) yg)*dx-Y[np][1];
			gaussy[np][i] = anorm*std::exp(-xx*xx/anorm2);
			gaussy_dip[np][i] = anormdip*std::exp(-xx*xx/anormdip2);
			grad_gaussy_dip[np][i] = -xx / sigmadipsq;

			zg = zc - ngdh + (i+1);
			indz[np][i] = zg - NPTS_Z * ((int) std::floor( ((double) zg) / ((double) NPTS_Z)));
			xx = ((double) zg)*dx-Y[np][2];
			gaussz[np][i] = anorm*std::exp(-xx*xx/anorm2);
			gaussz_dip[np][i] = anormdip*std::exp(-xx*xx/anormdip2);
			grad_gaussz_dip[np][i] = -xx / sigmadipsq;
		}
	}
	return;
}


void mono_and_dipole_force_distribution(
	double ***fx, double ***fy, double ***fz,
	double **GS, double **GA, double **F,
	double **gaussx, double **gaussy, double **gaussz,
	double **gaussx_dip, double **gaussy_dip, double **gaussz_dip,
	double **grad_gaussx_dip, double **grad_gaussy_dip,
	double **grad_gaussz_dip,
	int **indx, int **indy, int **indz,
	int local_x_start, int local_nx){
	//double *fx2, double *fy2, double *fz2){
	int np, i, j, k, ii, jj, kk;
	double xx, yy, zz, temp, temp2;
	double g11, g22, g33, g12, g21, g13, g31, g23, g32;

	for(np = 0; np < Np; np++) {
		g11 = GS[np][0] + GA[np][0];
		g22 = GS[np][1] + GA[np][1];
		g33 = GS[np][2] + GA[np][2];
		g12 = GS[np][3] + GA[np][3];
		g21 = GS[np][3] - GA[np][3];
		g13 = GS[np][4] + GA[np][4];
		g31 = GS[np][4] - GA[np][4];
		g23 = GS[np][5] + GA[np][5];
		g32 = GS[np][5] - GA[np][5];
		for(i = 0; i < ngd; i++) {
			ii = indx[np][i];
			if( (local_x_start <= ii) && ii < (local_x_start + local_nx)) {
				ii = ii-local_x_start;
				xx = grad_gaussx_dip[np][i];
				for(j = 0; j < ngd; j++) {
					jj = indy[np][j];
					yy = grad_gaussy_dip[np][j];
					for(k = 0; k < ngd; k++) {
						kk = indz[np][k];
						zz = grad_gaussz_dip[np][k];
						temp = gaussx[np][i] * gaussy[np][j] * gaussz[np][k];
						temp2 = gaussx_dip[np][i] * gaussy_dip[np][j]
						        * gaussz_dip[np][k];
						fx[ii][jj][kk] = fx[ii][jj][kk] + F[np][0]*temp
						                 + (g11*xx + g12*yy + g13*zz)*temp2;
						fy[ii][jj][kk] = fy[ii][jj][kk] + F[np][1]*temp
						                 + (g21*xx + g22*yy + g23*zz)*temp2;
						fz[ii][jj][kk] = fz[ii][jj][kk] + F[np][2]*temp
						                 + (g31*xx + g32*yy + g33*zz)*temp2;
					}
				}
			}
		}
	}
	return;
}

void flow_solve(double ***fx, double ***fy, double ***fz,
                double ***ux, double ***uy, double ***uz,
                double *qx, double *qy, double *qpad,
                double *qxsq, double *qysq, double *qpadsq,
                int local_nyt, int local_y_start, int pad){

	double totalpts, smallx = 1e-18;
	double q1, q2, q3, qq, f1, f2,f3;
	double *f1_p, *f2_p, *f3_p, *u1_p, *u2_p, *u3_p;
	double norm, kdotf;

	if(local_y_start == 0) {
		fx[0][0][0] = 0.0;
		fx[0][0][1] = 0.0;
		fy[0][0][0] = 0.0;
		fy[0][0][1] = 0.0;
		fz[0][0][0] = 0.0;
		fz[0][0][1] = 0.0;
	}

	//cout << "XB" << endl;

	//MPI_Barrier(MPI_COMM_WORLD); // remove later

	totalpts = (double) (NPTS_X*NPTS_Y*NPTS_Z); // XXX try to openmp this??
	// Point to the beginning of the fx, fy, fz arrays
	f1_p = &fx[0][0][0];
	f2_p = &fy[0][0][0];
	f3_p = &fz[0][0][0];
	u1_p = &ux[0][0][0];
	u2_p = &uy[0][0][0];
	u3_p = &uz[0][0][0];
	for(int j = 0; j < local_nyt; j++) {
		//cout << "[" << j+local_y_start;
		q2 = qy[j+local_y_start];
		//cout  << "]" << endl;
		for(int i = 0; i < NPTS_X; i++) {
			//cout << "{" << j << "," << i << "}";
			q1 = qx[i];
			qq = qysq[j+local_y_start] + qxsq[i] + smallx;
			for(int k = 0; k < pad; k++) {
				//cout << "{" << j << "," << i << "," << k << "}" <<endl;
				q3 = qpad[k];
				//cout << "." <<endl;
				f1 = *f1_p++;
				f2 = *f2_p++;
				f3 = *f3_p++;
				// f1 = value at memory address f1_p, then advance address
				// f1_p by one. Same for f2, f3.
				// This is equivalent to saying
				//   f1 = *(&fx[0][0][0] + j*NPTS_X*pad + i*pad + k);
				//   f2 = *(&fy[0][0][0] + j*NPTS_X*pad + i*pad + k);
				//   f3 = *(&fz[0][0][0] + j*NPTS_X*pad + i*pad + k);
				// which in turn is equivalent to fx[j][i][k],
				// I think, but allows overflowing.
				// Compare with FFTW docs to be convinced, or uncomment:
				// cout << "(" << fx[j][i][k] << "==" << f1 << "?)" << endl;

				/*
				   if (i<10 && j<10 && k<10){
				    // These won't agree when j!=0 and that's OK I think.
				    cout << "(" << j << "," << i << "," << k << " " << fx[j][i][k] << "==" << f1 << "?)" << endl;
				   }
				 */
				//cout << "." <<endl;
				norm = 1.0/(qq + qpadsq[k]);
				kdotf = (q1*f1+q2*f2+q3*f3)*norm;
				//cout << "." <<endl;
				// Equivalent to ux[j][i][k] as above.
				// Note i++ = 1 is equivalent to i=1; i++.
				*u1_p++ = norm*(f1-q1*(kdotf))/totalpts;
				*u2_p++ = norm*(f2-q2*(kdotf))/totalpts;
				*u3_p++ = norm*(f3-q3*(kdotf))/totalpts;
			}
		}
	}

	//cout << "XC" << endl;

	if(local_y_start == 0) {
		ux[0][0][0] = 0.0;
		ux[0][0][1] = 0.0;
		uy[0][0][0] = 0.0;
		uy[0][0][1] = 0.0;
		uz[0][0][0] = 0.0;
		uz[0][0][1] = 0.0;
	}
	return;
}



void particle_velocities(double ***ux, double ***uy, double ***uz, double **VTEMP,
                         double **gaussx, double **gaussy, double **gaussz,
                         int **indx, int **indy, int **indz,
                         int local_x_start, int local_nx, double dx){

	int np, i, j, k, ii, jj, kk;
	double norm, temp, x;

	norm = dx*dx*dx;

	for(np = 0; np < Np; np++) {
		for(i = 0; i < ngd; i++) {
			ii = indx[np][i];
			x = dx * ((double) ii);
			if((local_x_start <= ii) && ii < (local_x_start + local_nx)) {
				ii = ii-local_x_start;
				for(j = 0; j < ngd; j++) {
					jj = indy[np][j];
					for(k = 0; k < ngd; k++) {
						kk = indz[np][k];
						temp = gaussx[np][i]*gaussy[np][j]*gaussz[np][k]*norm;
						VTEMP[np][0] = VTEMP[np][0]+ux[ii][jj][kk]*temp;
						VTEMP[np][1] = VTEMP[np][1]+uy[ii][jj][kk]*temp;
						VTEMP[np][2] = VTEMP[np][2]+uz[ii][jj][kk]*temp;
					}
				}
			}
		}
	}
	return;
}



void particle_rotations(double ***ux, double ***uy, double ***uz, double **WTEMP,
                        double **gaussx_dip, double **gaussy_dip, double **gaussz_dip,
                        double **grad_gaussx_dip, double **grad_gaussy_dip, double **grad_gaussz_dip,
                        int **indx, int **indy, int **indz,
                        int local_x_start, int local_nx, double dx){

	int np, i, j, k, ii, jj, kk;
	double x, xx, yy, zz, temp, norm;

	norm = dx*dx*dx;

	for(np = 0; np < Np; np++) {
		for(i = 0; i < ngd; i++) {
			ii = indx[np][i];
			x = dx * ((double) ii);
			if((local_x_start <= ii) && ii < (local_x_start + local_nx)) {
				ii = ii-local_x_start;
				xx = grad_gaussx_dip[np][i];
				for(j = 0; j < ngd; j++) {
					jj = indy[np][j];
					yy = grad_gaussy_dip[np][j];
					for(k = 0; k < ngd; k++) {
						kk = indz[np][k];
						zz = grad_gaussz_dip[np][k];
						temp = gaussx_dip[np][i]*gaussy_dip[np][j]*gaussz_dip[np][k]*norm;
						WTEMP[np][0] = WTEMP[np][0]-0.5*(uz[ii][jj][kk]*yy - uy[ii][jj][kk]*zz)*temp;
						WTEMP[np][1] = WTEMP[np][1]-0.5*(ux[ii][jj][kk]*zz - uz[ii][jj][kk]*xx)*temp;
						WTEMP[np][2] = WTEMP[np][2]-0.5*(uy[ii][jj][kk]*xx - ux[ii][jj][kk]*yy)*temp;
					}
				}
			}
		}
	}
	return;
}


void FCMStokeslet(double **Y, double **F,double **T,double **V,double **W){
	double a1 = 1./(6*PI*a*mu);
	double a2 = 1./(8*PI*a*a*a*mu);
	for(int np = 0; np < Np; ++np) {

		V[np][0] = a1*F[np][0];
		V[np][1] = a1*F[np][1];
		V[np][2] = a1*F[np][2];

		W[np][0] = a2*T[np][0];
		W[np][1] = a2*T[np][1];
		W[np][2] = a2*T[np][2];
	}

	//add stokeslet
	for(int i = 0; i < Np; ++i) {
		vec Xi(3); // position
		vec Xj(3); // position
		vec Xij(3); // displacement
		vec Fi(3); // Force
		vec Fj(3); // Force
		vec Vi(3); // Velocity
		vec Vj(3); // Velocity
		vec Rhat(3); // direction
		double r; // distance
		double RhatDotTmp; // something projected onto displacement

		Xi <<  Y[i][0] <<  Y[i][1] << Y[i][2];
		Fi <<  F[i][0] <<  F[i][1] << F[i][2];
		Vi <<  V[i][0] <<  V[i][1] << V[i][2];

		for(int j = 0; j < Np; ++j) {
			if(i != j) {

				Xj <<  Y[j][0] <<  Y[j][1] << Y[j][2];
				Fj <<  F[j][0] <<  F[j][1] << F[j][2];
				Vj <<  V[j][0] <<  V[j][1] << V[j][2];

				Xij  = (Xi-Xj);
				r    = norm(Xij);
				Rhat = Xij/r;

				double a3 = 1./(8.*PI*r);
				RhatDotTmp = dot(Fj,Rhat);
				Vi += a3*Fj + a3*RhatDotTmp*Rhat;
			}
		}
		V[i][0] = Vi(0);
		V[i][1] = Vi(1);
		V[i][2] = Vi(2);
	}
	return;
}

void Fcm::friction_mobility_solve(){
	double a1 = 1./(6*PI*a*mu);
	double a2 = 1./(8*PI*a*a*a*mu);
	for(int np = 0; np < Np; ++np) {

		this->V[np][0] = a1*this->F[np][0];
		this->V[np][1] = a1*this->F[np][1];
		this->V[np][2] = a1*this->F[np][2];

		this->W[np][0] = a2*this->T[np][0];
		this->W[np][1] = a2*this->T[np][1];
		this->W[np][2] = a2*this->T[np][2];
	}
	MPI_Barrier(MPI_COMM_WORLD); // remove later. yes actually remove it later! 14/1/20
	return;
}

void Fcm::assign_filament_data(const std::vector<Filament>& filaments){
	// only run by rank 0

	// copies and scales torque by LfcmBox
	for(int j = 0; j<Np; ++j) {

		// search filament index and bead index within filament ==============
		int n = 0; // index of filament
		int lengthbefore = 0; // index of first element in arrays

		while (lengthbefore<=j && n<filaments.size()) {
			lengthbefore += filaments[n].length();
			++n;
		}
		--n; // filament index
		lengthbefore -= filaments[n].length();
		int i = j - lengthbefore; // swimmer index

		//cout << "Particle " << j << " is Fil " << n << " Bead " << i << endl;

		this->Y[j][0] = std::fmod(filaments[n].beads[i].Xs(0),LfcmBox_x);
		this->Y[j][1] = std::fmod(filaments[n].beads[i].Xs(1),LfcmBox_y);
		this->Y[j][2] = std::fmod(filaments[n].beads[i].Xs(2),LfcmBox_z);

		// forces
		this->F[j][0] = filaments[n].beads[i].F(0);
		this->F[j][1] = filaments[n].beads[i].F(1);
		this->F[j][2] = filaments[n].beads[i].F(2);

		// torques
		this->T[j][0] = filaments[n].beads[i].TAU(0);
		this->T[j][1] = filaments[n].beads[i].TAU(1);
		this->T[j][2] = filaments[n].beads[i].TAU(2);
	}

	return;
}

void Fcm::mpi_broadcast(){
	MPI_Bcast(&this->Y[0][0], 3*Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&this->F[0][0], 3*Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&this->T[0][0], 3*Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
