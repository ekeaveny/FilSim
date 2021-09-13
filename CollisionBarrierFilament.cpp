#include "multi_filament_header.hpp"
#include "CollisionBarrierFilament.hpp"
#include "c_array_functions.hpp"

// inclusion guards
//#ifndef FILCOLBAR_IS_INCLUDED
//#define FILCOLBAR_IS_INCLUDED

// uses the filament structure
// sfs14@ic.ac.uk. started 25/09/18

using namespace arma;

// Collision barrier repulsion force between beads. Works for periodic box and
// does not include repulsion between neighbouring beads in filament.
// Note 1: strength F^S must be greater than largest K_B/L^2.
// Note 2: This is called by each process.
/*void collision_barrier_linkedlist(std::vector<Filament>& filaments,
                                  int *map, int *head, int *list,
                                  int loc_cell_start, int loc_cell_stop,
                                  int *global_segment_to_local_filament_number,
                                  int *global_segment_to_local_segment_number){*/
//void linked_list::apply_collision_barrier(std::vector<Filament>& filaments){
void LinkedList::apply_collision_barrier(Fcm& fcm, int nt){
	int i, j, jcello, jcell;
	int i_fil, i_seg, j_fil, j_seg;
	double xi, yi, zi, fxi, fyi, fzi, xij, yij, zij, rij, rijsq, FXij, FYij, FZij;
	double temp, temp2;
	double LfcmBoxhalf_x = LfcmBox_x*.5;
	double LfcmBoxhalf_y = LfcmBox_y*.5;
	double LfcmBoxhalf_z = LfcmBox_z*.5;
	const double Rnm  = 2*a;
	const double Rnm2 = Rnm*Rnm;
	const double chiRnm2 = 1.21*Rnm2;
	double prefac;
	// Collision barrier strength F^S as in paper, (21), divided by 2a.
	// Must have F^S > K_B/L^2.
	double FS_over_2a = Barrier_FS_over_2a; // from config.
	double force_cap_distance = Barrier_distance_cap; // from config //sqrt(chiRnm2); //
	double force_cap_distance2 = force_cap_distance*force_cap_distance; //sqrt(chiRnm2); //
	double min_dist2 = 100000;
	double min_dist2_over_processors = 100000;
	//double min_dist_force_x = 0.0;
	//double min_dist_force_y = 0.0;
	//double min_dist_force_z = 0.0;
	//bool min_force = false;

	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	for (int i=0; i < Np; ++i) {
		for (int j=0; j < 3; ++j) {
			this->F_on_node[i][j] = 0;
			this->F_total[i][j] = 0;
		}
	}

	for(int icell = this->loc_cell_start; icell < this->loc_cell_stop; icell++) {
		i = this->head[icell];
		while(i > -1) {
			i_fil = this->global_segment_to_local_filament_number[i];
			i_seg = this->global_segment_to_local_segment_number[i];
			//cout << myrank << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]" << endl;
			xi = fcm.Y[i][0];
			yi = fcm.Y[i][1];
			zi = fcm.Y[i][2];
			fxi = this->F_on_node[i][0];
			fyi = this->F_on_node[i][1];
			fzi = this->F_on_node[i][2];
			j = this->list[i];
			while(j > -1) {
				j_fil = this->global_segment_to_local_filament_number[j];
				j_seg = this->global_segment_to_local_segment_number[j];
				//cout << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]" << ", j = " << j << " = [" << j_fil << "," << j_seg << "]" << endl;
				//cin.get();

				// No barrier between adjacent segments in same filament
				// unless we are using Stokes drag ('friction') mobility solve
				// (because nothing to stop them overlapping otherwise)
				if (!(j_fil == i_fil && abs(j_seg-i_seg)==1) || nt < initialFrictionSteps) {
					xij = xi - fcm.Y[j][0];
					yij = yi - fcm.Y[j][1];
					zij = zi - fcm.Y[j][2];

					// forces -L/2 < x < L/2 where L=Lfcmbox
					//xij = xij - LfcmBox_x * ((double) ((int) (xij/LfcmBoxhalf_x)));
					//yij = yij - LfcmBox_y * ((double) ((int) (yij/LfcmBoxhalf_y)));
					//zij = zij - LfcmBox_z * ((double) ((int) (zij/LfcmBoxhalf_z)));
					xij = xij - LfcmBox_x * floor(xij/LfcmBox_x + 0.5);
					yij = yij - LfcmBox_y * floor(yij/LfcmBox_y + 0.5);
					zij = zij - LfcmBox_z * floor(zij/LfcmBox_z + 0.5);

					rijsq=xij*xij+yij*yij+zij*zij;
					/*
					   if (((i==0 && j==22)||(i==22&&j==0))){
					    cout << "0_22 dist is " << std::sqrt(rijsq) << endl;
					   }
					 */

					//min_force = false;
					if (rijsq < min_dist2) {
						min_dist2 = rijsq;
						//min_dist_force_x = 0;
						//min_dist_force_y = 0;
						//min_dist_force_z = 0;
						//min_force = true;
					}

					if(rijsq < chiRnm2) {

						if (rijsq < force_cap_distance2) {
							// Cap magnitude at when rij = 1
							rij = sqrt(rijsq);
							xij = xij/rij*force_cap_distance;
							yij = yij/rij*force_cap_distance;
							zij = zij/rij*force_cap_distance;
							rijsq = force_cap_distance2;
						}

						prefac  =  chiRnm2 - rijsq;
						prefac /= (chiRnm2 - Rnm2);
						prefac *= prefac;
						prefac *= prefac; // fourth power!
						prefac *= FS_over_2a;

						FXij = prefac*xij;
						FYij = prefac*yij;
						FZij = prefac*zij;

						fxi += FXij;
						fyi += FYij;
						fzi += FZij;
						this->F_on_node[j][0] -= FXij;
						this->F_on_node[j][1] -= FYij;
						this->F_on_node[j][2] -= FZij;

						//if (min_force) {
						//	min_dist_force_x = FXij;
						//	min_dist_force_y = FYij;
						//	min_dist_force_z = FZij;
						//}

						//if(myrank == 0){
						//if (i_fil == 10 && j_fil == 7){
						//cout << "Colbar force " << i_fil << "," << i_seg << "-" << j_fil << "," << j_seg << " (dist " << sqrt(rijsq) << ")" << " force+ on " << i_fil << ": " << FXij << "," << FYij << "," << FZij << endl;
						//cout << "Collision barrier between [" << i_fil << "][" << i_seg << "] and [" << j_fil << "][" << j_seg << "], dist2=" << rijsq << endl;
						//}
						//}
					}
				}
				j = this->list[j];
			}

			jcello = 13*icell;
			for(int nabor = 0; nabor < 13; nabor++) {
				jcell = this->map[jcello + nabor];
				j = this->head[jcell];
				while(j > -1) {
					j_fil = this->global_segment_to_local_filament_number[j];
					j_seg = this->global_segment_to_local_segment_number[j];
					//cout << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]," << "Cell " << jcell << ": j = " << j << " = [" << j_fil << "," << j_seg << "]" << endl;

					// No barrier between adjacent segments in same filament
					// unless we are using Stokes drag ('friction') mobility solve
					// (because nothing to stop them overlapping otherwise)
					if (!(j_fil == i_fil && abs(j_seg-i_seg)==1) || nt < initialFrictionSteps) {
						xij = xi - fcm.Y[j][0];
						yij = yi - fcm.Y[j][1];
						zij = zi - fcm.Y[j][2];

						// forces -L/2 < x < L/2 where L=Lfcmbox
						//xij = xij - LfcmBox_x * ((double) ((int) (xij/LfcmBoxhalf_x)));
						//yij = yij - LfcmBox_y * ((double) ((int) (yij/LfcmBoxhalf_y)));
						//zij = zij - LfcmBox_z * ((double) ((int) (zij/LfcmBoxhalf_z)));
						xij = xij - LfcmBox_x * floor(xij/LfcmBox_x + 0.5);
						yij = yij - LfcmBox_y * floor(yij/LfcmBox_y + 0.5);
						zij = zij - LfcmBox_z * floor(zij/LfcmBox_z + 0.5);

						rijsq=xij*xij+yij*yij+zij*zij;

						/*
						   if (((i==0 && j==22)||(i==22&&j==0))){
						    cout << i << "-" << j << " dist is " << std::sqrt(rijsq)
						        << "|" << fcm.Y[i][0] << " - " << fcm.Y[j][0] << " = " << xij << "? "
						        << "|" << fcm.Y[i][1] << " - " << fcm.Y[j][1] << " = " << yij << "? "
						        << "|" << fcm.Y[i][2] << " - " << fcm.Y[j][2] << " = " << zij << "? "
						        << endl;
						   }
						 */
						//min_force = false;
						if (rijsq < min_dist2) {
							min_dist2 = rijsq;
							//min_dist_force_x = 0;
							//min_dist_force_y = 0;
							//min_dist_force_z = 0;
							//min_force = true;
						}

						if(rijsq < chiRnm2) {

							if (rijsq < force_cap_distance2) {
								// Cap magnitude at when rij = force_cap_distance
								rij = sqrt(rijsq);
								xij = xij/rij*force_cap_distance;
								yij = yij/rij*force_cap_distance;
								zij = zij/rij*force_cap_distance;
								rijsq = force_cap_distance2;
							}

							prefac  =  chiRnm2 - rijsq;
							prefac /= (chiRnm2 - Rnm2);
							prefac *= prefac;
							prefac *= prefac; // fourth power!
							prefac *= FS_over_2a;

							FXij = prefac*xij;
							FYij = prefac*yij;
							FZij = prefac*zij;

							fxi += FXij;
							fyi += FYij;
							fzi += FZij;

							this->F_on_node[j][0] -= FXij;
							this->F_on_node[j][1] -= FYij;
							this->F_on_node[j][2] -= FZij;

							//if (min_force) {
							//	min_dist_force_x = FXij;
							//	min_dist_force_y = FYij;
							//	min_dist_force_z = FZij;
							//}
							//if(myrank == 0){
							//if (i_fil == 10 && j_fil == 7){
							//cout << "Colbar forc " << i_fil << "," << i_seg << "-" << j_fil << "," << j_seg << " (dist " << sqrt(rijsq) << ")" << " force+ on " << i_fil << ": " << FXij << "," << FYij << "," << FZij << endl;
							//cout << "Collision barrier between [" << i_fil << "][" << i_seg << "] and [" << j_fil << "][" << j_seg << "], dist2=" << rijsq << endl;
							//}							//}
						}
					}
					j = this->list[j];
				}
			}
			this->F_on_node[i][0] = fxi;
			this->F_on_node[i][1] = fyi;
			this->F_on_node[i][2] = fzi;
			i = this->list[i];
		}
	}

	MPI_Allreduce(&min_dist2, &this->min_dist2, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	/*if (abs(this->min_dist2 - min_dist2_over_processors)<1e-6){
	    cout << "Min_dist: " << sqrt(min_dist2) << endl; // << ", force = [" << min_dist_force_x << ", " << min_dist_force_y << ", " << min_dist_force_z << "]" << endl;
	    this->min_dist2 = min_dist2;
	   }*/

	MPI_Allreduce(&this->F_on_node[0][0], &this->F_total[0][0], 3*Np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//if (myrank == 0){
	for (int i=0; i<Np; ++i) {
		for (int j=0; j<3; ++j) {
			fcm.F[i][j] += this->F_total[i][j];
			//if (abs(this->F_total[i][j]) > 1e-10){
			//	cout << "F_total["<<i<<"]["<<j<<"] = " << this->F_total[i][j] << " ==> F["<<i<<"]["<<j<<"] = " << fcm.F[i][j] << endl;
			//}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//}

	return;
}


void LinkedList::initialise(int myrank, int totalnodes,
                            std::vector<Filament>& filaments){

	Rref = 2.5*a; // Ideal cell width. In practice will be slighly larger ...
	M_X = (int) (LfcmBox_x/Rref); // ... because of the (int) here.
	M_Y = (int) (LfcmBox_y/Rref); // ... because of the (int) here.
	M_Z = (int) (LfcmBox_z/Rref); // ... because of the (int) here.
	ncell = M_X*M_Y*M_Z;
	mapsize=13*ncell;

	loc_cell_start = (myrank) * ((int) (((double) ncell) / ((double) totalnodes)));
	loc_cell_stop = (myrank+1) * ((int) (((double) ncell) / ((double) totalnodes)));

	if(myrank == totalnodes-1) {
		loc_cell_stop = ncell;
	}

	map = int_array_1D(mapsize);
	// One head and list to keep track of filament number, another to keep
	// track of segment (bead) number.
	head = int_array_1D(ncell);
	list = int_array_1D(Np);
	bulkmap(map, M_X, M_Y, M_Z);

	// Set up global segment number <-> [fil,seg] arrays.
	set_global_segment_number(filaments); /*, global_segment_number,
	                                         global_segment_to_local_filament_number,
	                                         global_segment_to_local_segment_number); */

	Fref = 2.0/a;
}


int LinkedList::icell(int M_X, int M_Y, int M_Z, int x, int y, int z){
	double beta;
	int q;
	double cellix = (double) M_X;
	double celliy = (double) M_Y;
	beta=fmod((x+M_X), M_X)+fmod((y+M_Y), M_Y)*cellix+fmod((z+M_Z), M_Z)*cellix*celliy;
	q=(int)(beta);
	return q;
}


void LinkedList::bulkmap(int *map, int M_X, int M_Y, int M_Z){
	int imap=0, tempmap=0;
	int iz = 0, iy = 0, ix = 0;
	for(iz = 0; iz < M_Z; iz++) {
		for(iy = 0; iy < M_Y; iy++) {
			for(ix = 0; ix < M_X; ix++) {
				tempmap=icell(M_X, M_Y, M_Z, ix, iy, iz);
				imap=tempmap*13;
				map[imap]=icell(M_X, M_Y, M_Z, ix+1, iy, iz);
				map[imap+1]=icell(M_X, M_Y, M_Z, ix+1, iy+1, iz);
				map[imap+2]=icell(M_X, M_Y, M_Z, ix, iy+1, iz);
				map[imap+3]=icell(M_X, M_Y, M_Z, ix-1, iy+1, iz);
				map[imap+4]=icell(M_X, M_Y, M_Z, ix+1, iy, iz-1);
				map[imap+5]=icell(M_X, M_Y, M_Z, ix+1, iy+1, iz-1);
				map[imap+6]=icell(M_X, M_Y, M_Z, ix, iy+1, iz-1);
				map[imap+7]=icell(M_X, M_Y, M_Z, ix-1, iy+1, iz-1);
				map[imap+8]=icell(M_X, M_Y, M_Z, ix+1, iy, iz+1);
				map[imap+9]=icell(M_X, M_Y, M_Z, ix+1, iy+1, iz+1);
				map[imap+10]=icell(M_X, M_Y, M_Z, ix, iy+1, iz+1);
				map[imap+11]=icell(M_X, M_Y, M_Z, ix-1, iy+1, iz+1);
				map[imap+12]=icell(M_X, M_Y, M_Z, ix, iy, iz+1);
			}
		}
	}
	return;
}

void LinkedList::set_global_segment_number(vector<Filament>& filaments){
	// Run by every rank but only rank 0 has the correct filament number and
	// length data, hence the broadcast at the end from rank 0 to all ranks.

	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int g = 0;
	for(int i_fil=0; i_fil < filaments.size(); ++i_fil) {
		//cout << myrank << " " << i_fil << " has length " << filaments[i_fil].myNworm << endl;
		for(int i_seg=0; i_seg < filaments[i_fil].length(); ++i_seg) {
			this->global_segment_number[i_fil][i_seg] = g;
			this->global_segment_to_local_filament_number[g] = i_fil;
			this->global_segment_to_local_segment_number[g] = i_seg;
			++g;
		}
	}

	MPI_Bcast(&this->global_segment_number[0][0], Nsw*max_Nworm, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&this->global_segment_to_local_filament_number[0], Np, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&this->global_segment_to_local_segment_number[0], Np, MPI_INT, 0, MPI_COMM_WORLD);
	return;
}


void LinkedList::link(vector<Filament>& filaments){
	// Only run by rank 0
	int index;
	double xr, yr, zr;
	double celli, cellj, cellk;
	for(int i = 0; i<this->ncell; i++) {
		this->head[i] = -1;
	}
	celli = (double) (this->M_X);
	cellj = (double) (this->M_Y);
	cellk = (double) (this->M_Z);
	int global_segment_number = 0;

	int large_number = 1e5;
	for(int i_fil=0; i_fil < filaments.size(); ++i_fil) {
		for(int i_seg=0; i_seg < filaments[i_fil].myNworm; ++i_seg) {
			//cout << "(" << endl;
			xr = filaments[i_fil].beads[i_seg].Xs(0);
			yr = filaments[i_fil].beads[i_seg].Xs(1);
			zr = filaments[i_fil].beads[i_seg].Xs(2);

			// BEGIN ADDED 28/6/19
			// forces 0 < x < L where L=Lfcmbox
			xr = xr - LfcmBox_x * ((double) ((int) (xr/LfcmBox_x + large_number) - large_number));
			yr = yr - LfcmBox_y * ((double) ((int) (yr/LfcmBox_y + large_number) - large_number));
			zr = zr - LfcmBox_z * ((double) ((int) (zr/LfcmBox_z + large_number) - large_number));
			// END ADDED 28/6/19

			// Check for explosion
			if (abs(xr) > large_number || abs(yr) > large_number || abs(zr) > large_number) {
				cout << endl;
				cout << "ERROR: Coordinates too large to place in periodic box. ";
				cout << "Normally due to explosion. Segfault imminent." << endl;
				cout << endl;
			}

			//cout << "I" << i_fil << "," << i_seg << "," << endl;
			index = (int)((xr/LfcmBox_x)*celli)
			        + (int)((yr/LfcmBox_y)*cellj)*M_X
			        + (int)((zr/LfcmBox_z)*cellk)*M_X*M_Y;
			//cout << "Il" << xr << "," << LfcmBox_x << endl;
			//cout << "Im" << (xr/LfcmBox_x) << "," << celli << "," << M_X << "," << M_Y << "," << (yr/LfcmBox_y) << "," << (zr/LfcmBox_z) << endl;
			//cout << "J" << index << "," << global_segment_number << "," << endl;
			//cout << "K" << list[global_segment_number] << "," << endl;
			//cout << "L" << head[index] << "" << endl;
			this->list[global_segment_number]=this->head[index];
			this->head[index]=global_segment_number;
			++global_segment_number;
			//cout << ")" << endl;
		}
	}
	/*
	   int l = 0;
	   for (int k=0; k<M*M*M; ++k){
	    if (this->head[k] != -1){
	        cout << "(" << l << ") head[" << k << "] = " << this->head[k];
	        cout << ". Next on list = " << this->list[this->head[k]];
	        cout << endl;
	 ++l;
	    }
	   }
	 */
	return;
}

void LinkedList::mpi_broadcast(){
	MPI_Bcast(&this->head[0], this->ncell, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&this->list[0], Np, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&this->map[0], this->mapsize, MPI_INT, 0, MPI_COMM_WORLD);
}

void LinkedList::apply_collision_barrier_not_swimmer(Fcm& fcm, int nt){
	int i, j, jcello, jcell;
	int i_fil, i_seg, j_fil, j_seg;
	double xi, yi, zi, fxi, fyi, fzi, xij, yij, zij, rij, rijsq, FXij, FYij, FZij;
	double temp, temp2;
	double LfcmBoxhalf_x = LfcmBox_x*.5;
	double LfcmBoxhalf_y = LfcmBox_y*.5;
	double LfcmBoxhalf_z = LfcmBox_z*.5;
	const double Rnm  = 2*a;
	const double Rnm2 = Rnm*Rnm;
	const double chiRnm2 = 1.21*Rnm2;
	double prefac;
	// Collision barrier strength F^S as in paper, (21), divided by 2a.
	// Must have F^S > K_B/L^2.
	double FS_over_2a = Barrier_FS_over_2a; // from config.
	double force_cap_distance = Barrier_distance_cap; // from config //sqrt(chiRnm2); //
	double force_cap_distance2 = force_cap_distance*force_cap_distance; //sqrt(chiRnm2); //
	double min_dist2 = 100000;
	double min_dist2_over_processors = 100000;
	//double min_dist_force_x = 0.0;
	//double min_dist_force_y = 0.0;
	//double min_dist_force_z = 0.0;
	//bool min_force = false;

	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	for (int i=0; i < Np; ++i) {
		for (int j=0; j < 3; ++j) {
			this->F_on_node[i][j] = 0;
			this->F_total[i][j] = 0;
		}
	}

	for(int icell = this->loc_cell_start; icell < this->loc_cell_stop; icell++) {
		i = this->head[icell];
		while(i > -1) {
			i_fil = this->global_segment_to_local_filament_number[i];
			i_seg = this->global_segment_to_local_segment_number[i];
			//cout << myrank << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]" << endl;
			xi = fcm.Y[i][0];
			yi = fcm.Y[i][1];
			zi = fcm.Y[i][2];
			fxi = this->F_on_node[i][0];
			fyi = this->F_on_node[i][1];
			fzi = this->F_on_node[i][2];
			j = this->list[i];
			while(j > -1) {
				j_fil = this->global_segment_to_local_filament_number[j];
				j_seg = this->global_segment_to_local_segment_number[j];
				//cout << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]" << ", j = " << j << " = [" << j_fil << "," << j_seg << "]" << endl;
				//cin.get();

				// No barriers between swimmer & network
				if((i_fil == 0) == (j_fil == 0)) {

					// No barrier between adjacent segments in same filament
					// unless we are using Stokes drag ('friction') mobility solve
					// (because nothing to stop them overlapping otherwise)
					if (!(j_fil == i_fil && abs(j_seg-i_seg)==1) || nt < initialFrictionSteps) {
						xij = xi - fcm.Y[j][0];
						yij = yi - fcm.Y[j][1];
						zij = zi - fcm.Y[j][2];

						// forces -L/2 < x < L/2 where L=Lfcmbox
						//xij = xij - LfcmBox_x * ((double) ((int) (xij/LfcmBoxhalf_x)));
						//yij = yij - LfcmBox_y * ((double) ((int) (yij/LfcmBoxhalf_y)));
						//zij = zij - LfcmBox_z * ((double) ((int) (zij/LfcmBoxhalf_z)));
						xij = xij - LfcmBox_x * floor(xij/LfcmBox_x + 0.5);
						yij = yij - LfcmBox_y * floor(yij/LfcmBox_y + 0.5);
						zij = zij - LfcmBox_z * floor(zij/LfcmBox_z + 0.5);

						rijsq=xij*xij+yij*yij+zij*zij;
						/*
						   if (((i==0 && j==22)||(i==22&&j==0))){
						    cout << "0_22 dist is " << std::sqrt(rijsq) << endl;
						   }
						 */

						//min_force = false;
						if (rijsq < min_dist2) {
							min_dist2 = rijsq;
							//min_dist_force_x = 0;
							//min_dist_force_y = 0;
							//min_dist_force_z = 0;
							//min_force = true;
						}

						if(rijsq < chiRnm2) {

							if (rijsq < force_cap_distance2) {
								// Cap magnitude at when rij = 1
								rij = sqrt(rijsq);
								xij = xij/rij*force_cap_distance;
								yij = yij/rij*force_cap_distance;
								zij = zij/rij*force_cap_distance;
								rijsq = force_cap_distance2;
							}

							prefac  =  chiRnm2 - rijsq;
							prefac /= (chiRnm2 - Rnm2);
							prefac *= prefac;
							prefac *= prefac; // fourth power!
							prefac *= FS_over_2a;

							FXij = prefac*xij;
							FYij = prefac*yij;
							FZij = prefac*zij;

							fxi += FXij;
							fyi += FYij;
							fzi += FZij;
							this->F_on_node[j][0] -= FXij;
							this->F_on_node[j][1] -= FYij;
							this->F_on_node[j][2] -= FZij;

							//if (min_force) {
							//	min_dist_force_x = FXij;
							//	min_dist_force_y = FYij;
							//	min_dist_force_z = FZij;
							//}

							//if(myrank == 0){
							//if (i_fil == 10 && j_fil == 7){
							//cout << "Colbar force " << i_fil << "," << i_seg << "-" << j_fil << "," << j_seg << " (dist " << sqrt(rijsq) << ")" << " force+ on " << i_fil << ": " << FXij << "," << FYij << "," << FZij << endl;
							//cout << "Collision barrier between [" << i_fil << "][" << i_seg << "] and [" << j_fil << "][" << j_seg << "], dist2=" << rijsq << endl;
							//}
							//}
						}
					}
				}
				j = this->list[j];
			}

			jcello = 13*icell;
			for(int nabor = 0; nabor < 13; nabor++) {
				jcell = this->map[jcello + nabor];
				j = this->head[jcell];
				while(j > -1) {
					j_fil = this->global_segment_to_local_filament_number[j];
					j_seg = this->global_segment_to_local_segment_number[j];
					//cout << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]," << "Cell " << jcell << ": j = " << j << " = [" << j_fil << "," << j_seg << "]" << endl;

					// No barriers between swimmer & network (!= equals XOR here)
					if((i_fil == 0) == (j_fil == 0)) {

						// No barrier between adjacent segments in same filament
						// unless we are using Stokes drag ('friction') mobility solve
						// (because nothing to stop them overlapping otherwise)
						if (!(j_fil == i_fil && abs(j_seg-i_seg)==1) || nt < initialFrictionSteps) {
							xij = xi - fcm.Y[j][0];
							yij = yi - fcm.Y[j][1];
							zij = zi - fcm.Y[j][2];

							// forces -L/2 < x < L/2 where L=Lfcmbox
							//xij = xij - LfcmBox_x * ((double) ((int) (xij/LfcmBoxhalf_x)));
							//yij = yij - LfcmBox_y * ((double) ((int) (yij/LfcmBoxhalf_y)));
							//zij = zij - LfcmBox_z * ((double) ((int) (zij/LfcmBoxhalf_z)));
							xij = xij - LfcmBox_x * floor(xij/LfcmBox_x + 0.5);
							yij = yij - LfcmBox_y * floor(yij/LfcmBox_y + 0.5);
							zij = zij - LfcmBox_z * floor(zij/LfcmBox_z + 0.5);

							rijsq=xij*xij+yij*yij+zij*zij;

							/*
							   if (((i==0 && j==22)||(i==22&&j==0))){
							    cout << i << "-" << j << " dist is " << std::sqrt(rijsq)
							        << "|" << fcm.Y[i][0] << " - " << fcm.Y[j][0] << " = " << xij << "? "
							        << "|" << fcm.Y[i][1] << " - " << fcm.Y[j][1] << " = " << yij << "? "
							        << "|" << fcm.Y[i][2] << " - " << fcm.Y[j][2] << " = " << zij << "? "
							        << endl;
							   }
							 */
							//min_force = false;
							if (rijsq < min_dist2) {
								min_dist2 = rijsq;
								//min_dist_force_x = 0;
								//min_dist_force_y = 0;
								//min_dist_force_z = 0;
								//min_force = true;
							}

							if(rijsq < chiRnm2) {

								if (rijsq < force_cap_distance2) {
									// Cap magnitude at when rij = force_cap_distance
									rij = sqrt(rijsq);
									xij = xij/rij*force_cap_distance;
									yij = yij/rij*force_cap_distance;
									zij = zij/rij*force_cap_distance;
									rijsq = force_cap_distance2;
								}

								prefac  =  chiRnm2 - rijsq;
								prefac /= (chiRnm2 - Rnm2);
								prefac *= prefac;
								prefac *= prefac; // fourth power!
								prefac *= FS_over_2a;

								FXij = prefac*xij;
								FYij = prefac*yij;
								FZij = prefac*zij;

								fxi += FXij;
								fyi += FYij;
								fzi += FZij;

								this->F_on_node[j][0] -= FXij;
								this->F_on_node[j][1] -= FYij;
								this->F_on_node[j][2] -= FZij;

								//if (min_force) {
								//	min_dist_force_x = FXij;
								//	min_dist_force_y = FYij;
								//	min_dist_force_z = FZij;
								//}
								//if(myrank == 0){
								//if (i_fil == 10 && j_fil == 7){
								//cout << "Colbar forc " << i_fil << "," << i_seg << "-" << j_fil << "," << j_seg << " (dist " << sqrt(rijsq) << ")" << " force+ on " << i_fil << ": " << FXij << "," << FYij << "," << FZij << endl;
								//cout << "Collision barrier between [" << i_fil << "][" << i_seg << "] and [" << j_fil << "][" << j_seg << "], dist2=" << rijsq << endl;
								//}							//}
							}
						}
					}
					j = this->list[j];
				}
			}
			this->F_on_node[i][0] = fxi;
			this->F_on_node[i][1] = fyi;
			this->F_on_node[i][2] = fzi;
			i = this->list[i];
		}
	}

	MPI_Allreduce(&min_dist2, &this->min_dist2, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	/*if (abs(this->min_dist2 - min_dist2_over_processors)<1e-6){
	    cout << "Min_dist: " << sqrt(min_dist2) << endl; // << ", force = [" << min_dist_force_x << ", " << min_dist_force_y << ", " << min_dist_force_z << "]" << endl;
	    this->min_dist2 = min_dist2;
	   }*/

	MPI_Allreduce(&this->F_on_node[0][0], &this->F_total[0][0], 3*Np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//if (myrank == 0){
	for (int i=0; i<Np; ++i) {
		for (int j=0; j<3; ++j) {
			fcm.F[i][j] += this->F_total[i][j];
			//if (abs(this->F_total[i][j]) > 1e-10){
			//	cout << "F_total["<<i<<"]["<<j<<"] = " << this->F_total[i][j] << " ==> F["<<i<<"]["<<j<<"] = " << fcm.F[i][j] << endl;
			//}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//}

	return;
}

void LinkedList::apply_collision_barrier_not_test_particles(Fcm& fcm, int nt){
	int i, j, jcello, jcell;
	int i_fil, i_seg, j_fil, j_seg;
	double xi, yi, zi, fxi, fyi, fzi, xij, yij, zij, rij, rijsq, FXij, FYij, FZij;
	double temp, temp2;
	double LfcmBoxhalf_x = LfcmBox_x*.5;
	double LfcmBoxhalf_y = LfcmBox_y*.5;
	double LfcmBoxhalf_z = LfcmBox_z*.5;
	const double Rnm  = 2*a;
	const double Rnm2 = Rnm*Rnm;
	const double chiRnm2 = 1.21*Rnm2;
	double prefac;
	// Collision barrier strength F^S as in paper, (21), divided by 2a.
	// Must have F^S > K_B/L^2.
	double FS_over_2a = Barrier_FS_over_2a; // from config.
	double force_cap_distance = Barrier_distance_cap; // from config //sqrt(chiRnm2); //
	double force_cap_distance2 = force_cap_distance*force_cap_distance; //sqrt(chiRnm2); //
	double min_dist2 = 100000;
	double min_dist2_over_processors = 100000;
	//double min_dist_force_x = 0.0;
	//double min_dist_force_y = 0.0;
	//double min_dist_force_z = 0.0;
	//bool min_force = false;

	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	for (int i=0; i < Np; ++i) {
		for (int j=0; j < 3; ++j) {
			this->F_on_node[i][j] = 0;
			this->F_total[i][j] = 0;
		}
	}

	for(int icell = this->loc_cell_start; icell < this->loc_cell_stop; icell++) {
		i = this->head[icell];
		while(i > -1) {
			i_fil = this->global_segment_to_local_filament_number[i];
			i_seg = this->global_segment_to_local_segment_number[i];
			//cout << myrank << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]" << endl;
			xi = fcm.Y[i][0];
			yi = fcm.Y[i][1];
			zi = fcm.Y[i][2];
			fxi = this->F_on_node[i][0];
			fyi = this->F_on_node[i][1];
			fzi = this->F_on_node[i][2];
			j = this->list[i];
			while(j > -1) {
				j_fil = this->global_segment_to_local_filament_number[j];
				j_seg = this->global_segment_to_local_segment_number[j];
				//cout << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]" << ", j = " << j << " = [" << j_fil << "," << j_seg << "]" << endl;
				//cin.get();

				// No barrier between adjacent segments in same filament
				// unless we are using Stokes drag ('friction') mobility solve
				// (because nothing to stop them overlapping otherwise)
				if ((!(j_fil == i_fil && abs(j_seg-i_seg)==1) || nt < initialFrictionSteps) && (i < Np-2 && j < Np-2)) {
					xij = xi - fcm.Y[j][0];
					yij = yi - fcm.Y[j][1];
					zij = zi - fcm.Y[j][2];

					// forces -L/2 < x < L/2 where L=Lfcmbox
					//xij = xij - LfcmBox_x * ((double) ((int) (xij/LfcmBoxhalf_x)));
					//yij = yij - LfcmBox_y * ((double) ((int) (yij/LfcmBoxhalf_y)));
					//zij = zij - LfcmBox_z * ((double) ((int) (zij/LfcmBoxhalf_z)));
					xij = xij - LfcmBox_x * floor(xij/LfcmBox_x + 0.5);
					yij = yij - LfcmBox_y * floor(yij/LfcmBox_y + 0.5);
					zij = zij - LfcmBox_z * floor(zij/LfcmBox_z + 0.5);

					rijsq=xij*xij+yij*yij+zij*zij;
					/*
					   if (((i==0 && j==22)||(i==22&&j==0))){
					    cout << "0_22 dist is " << std::sqrt(rijsq) << endl;
					   }
					 */

					//min_force = false;
					if (rijsq < min_dist2) {
						min_dist2 = rijsq;
						//min_dist_force_x = 0;
						//min_dist_force_y = 0;
						//min_dist_force_z = 0;
						//min_force = true;
					}

					if(rijsq < chiRnm2) {

						if (rijsq < force_cap_distance2) {
							// Cap magnitude at when rij = 1
							rij = sqrt(rijsq);
							xij = xij/rij*force_cap_distance;
							yij = yij/rij*force_cap_distance;
							zij = zij/rij*force_cap_distance;
							rijsq = force_cap_distance2;
						}

						prefac  =  chiRnm2 - rijsq;
						prefac /= (chiRnm2 - Rnm2);
						prefac *= prefac;
						prefac *= prefac; // fourth power!
						prefac *= FS_over_2a;

						FXij = prefac*xij;
						FYij = prefac*yij;
						FZij = prefac*zij;

						fxi += FXij;
						fyi += FYij;
						fzi += FZij;
						this->F_on_node[j][0] -= FXij;
						this->F_on_node[j][1] -= FYij;
						this->F_on_node[j][2] -= FZij;

						//if (min_force) {
						//	min_dist_force_x = FXij;
						//	min_dist_force_y = FYij;
						//	min_dist_force_z = FZij;
						//}

						//if(myrank == 0){
						//if (i_fil == 10 && j_fil == 7){
						//cout << "Colbar force " << i_fil << "," << i_seg << "-" << j_fil << "," << j_seg << " (dist " << sqrt(rijsq) << ")" << " force+ on " << i_fil << ": " << FXij << "," << FYij << "," << FZij << endl;
						//cout << "Collision barrier between [" << i_fil << "][" << i_seg << "] and [" << j_fil << "][" << j_seg << "], dist2=" << rijsq << endl;
						//}
						//}
					}
				}
				j = this->list[j];
			}

			jcello = 13*icell;
			for(int nabor = 0; nabor < 13; nabor++) {
				jcell = this->map[jcello + nabor];
				j = this->head[jcell];
				while(j > -1) {
					j_fil = this->global_segment_to_local_filament_number[j];
					j_seg = this->global_segment_to_local_segment_number[j];
					//cout << "Cell " << icell << ": i = " << i << " = [" << i_fil << "," << i_seg << "]," << "Cell " << jcell << ": j = " << j << " = [" << j_fil << "," << j_seg << "]" << endl;

					// No barrier between adjacent segments in same filament
					// unless we are using Stokes drag ('friction') mobility solve
					// (because nothing to stop them overlapping otherwise)
					if ((!(j_fil == i_fil && abs(j_seg-i_seg)==1) || nt < initialFrictionSteps) && (i < Np-2 && j < Np-2)) {
						xij = xi - fcm.Y[j][0];
						yij = yi - fcm.Y[j][1];
						zij = zi - fcm.Y[j][2];

						// forces -L/2 < x < L/2 where L=Lfcmbox
						//xij = xij - LfcmBox_x * ((double) ((int) (xij/LfcmBoxhalf_x)));
						//yij = yij - LfcmBox_y * ((double) ((int) (yij/LfcmBoxhalf_y)));
						//zij = zij - LfcmBox_z * ((double) ((int) (zij/LfcmBoxhalf_z)));
						xij = xij - LfcmBox_x * floor(xij/LfcmBox_x + 0.5);
						yij = yij - LfcmBox_y * floor(yij/LfcmBox_y + 0.5);
						zij = zij - LfcmBox_z * floor(zij/LfcmBox_z + 0.5);

						rijsq=xij*xij+yij*yij+zij*zij;

						/*
						   if (((i==0 && j==22)||(i==22&&j==0))){
						    cout << i << "-" << j << " dist is " << std::sqrt(rijsq)
						        << "|" << fcm.Y[i][0] << " - " << fcm.Y[j][0] << " = " << xij << "? "
						        << "|" << fcm.Y[i][1] << " - " << fcm.Y[j][1] << " = " << yij << "? "
						        << "|" << fcm.Y[i][2] << " - " << fcm.Y[j][2] << " = " << zij << "? "
						        << endl;
						   }
						 */
						//min_force = false;
						if (rijsq < min_dist2) {
							min_dist2 = rijsq;
							//min_dist_force_x = 0;
							//min_dist_force_y = 0;
							//min_dist_force_z = 0;
							//min_force = true;
						}

						if(rijsq < chiRnm2) {

							if (rijsq < force_cap_distance2) {
								// Cap magnitude at when rij = force_cap_distance
								rij = sqrt(rijsq);
								xij = xij/rij*force_cap_distance;
								yij = yij/rij*force_cap_distance;
								zij = zij/rij*force_cap_distance;
								rijsq = force_cap_distance2;
							}

							prefac  =  chiRnm2 - rijsq;
							prefac /= (chiRnm2 - Rnm2);
							prefac *= prefac;
							prefac *= prefac; // fourth power!
							prefac *= FS_over_2a;

							FXij = prefac*xij;
							FYij = prefac*yij;
							FZij = prefac*zij;

							fxi += FXij;
							fyi += FYij;
							fzi += FZij;

							this->F_on_node[j][0] -= FXij;
							this->F_on_node[j][1] -= FYij;
							this->F_on_node[j][2] -= FZij;

							//if (min_force) {
							//	min_dist_force_x = FXij;
							//	min_dist_force_y = FYij;
							//	min_dist_force_z = FZij;
							//}
							//if(myrank == 0){
							//if (i_fil == 10 && j_fil == 7){
							//cout << "Colbar forc " << i_fil << "," << i_seg << "-" << j_fil << "," << j_seg << " (dist " << sqrt(rijsq) << ")" << " force+ on " << i_fil << ": " << FXij << "," << FYij << "," << FZij << endl;
							//cout << "Collision barrier between [" << i_fil << "][" << i_seg << "] and [" << j_fil << "][" << j_seg << "], dist2=" << rijsq << endl;
							//}							//}
						}
					}
					j = this->list[j];
				}
			}
			this->F_on_node[i][0] = fxi;
			this->F_on_node[i][1] = fyi;
			this->F_on_node[i][2] = fzi;
			i = this->list[i];
		}
	}

	MPI_Allreduce(&min_dist2, &this->min_dist2, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	/*if (abs(this->min_dist2 - min_dist2_over_processors)<1e-6){
	    cout << "Min_dist: " << sqrt(min_dist2) << endl; // << ", force = [" << min_dist_force_x << ", " << min_dist_force_y << ", " << min_dist_force_z << "]" << endl;
	    this->min_dist2 = min_dist2;
	   }*/

	MPI_Allreduce(&this->F_on_node[0][0], &this->F_total[0][0], 3*Np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//if (myrank == 0){
	for (int i=0; i<Np; ++i) {
		for (int j=0; j<3; ++j) {
			fcm.F[i][j] += this->F_total[i][j];
			//if (abs(this->F_total[i][j]) > 1e-10){
			//	cout << "F_total["<<i<<"]["<<j<<"] = " << this->F_total[i][j] << " ==> F["<<i<<"]["<<j<<"] = " << fcm.F[i][j] << endl;
			//}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//}

	return;
}


//#endif
