// filament class. a filament consists of beads and lagrange multipliers
// representing their connections.
// it is inextensibile, but elastic.

// sfs14@ic.ac.uk, started 08/02/2018

// inclusion guards
#ifndef FILAMENT_IS_INCLUDED
#define FILAMENT_IS_INCLUDED

#include "config.hpp"
#include "Bead.hpp"

#include "FCMfunctions.hpp"
#include "CollisionBarrierFilament.hpp"
#include "spring_link_functions.hpp"

using namespace std;
using namespace arma;

class Filament {
private:

// variables that may differ between filaments

// strain-twist
vec myStrainTwist;

// undulation frequency & phase
double AngFrequency;
double Phase;

// wave Number, swimming helical alpha and beta
double Wavenumber;
double Alpha;
double Beta;
double L_filament;
double L_filament_long;
double myKAP;
double myC;
double myRhead;

// filament knows how many array elements before it are occupied by other filaments
// in Forc/Torque and state arrays.
int forceLengthFilled;     // occupied force array elements before filament
int stateLengthFilled;     // occupied state array elements before filament
int myNworm;     // number of beads making up filament
int myNlam;      // number of Lagrange multipliers


// mat myJapprox;
mat myJapproxInv;

// Lagrange multipliers
mat lam;
mat lam1;
mat lam2;

// Filament consist of BeadNumber beads


public:

Bead* beads;

friend class filamentJsolver;

// let error check access private member variables
friend bool ErrorCheck(const std::vector<Filament>& filaments, vec & Error,
                       const vec VX, const vec VY, const vec VZ,
                       const vec WX, const vec WY, const vec WZ, const int nt);

friend bool ErrorCheckFCM(const std::vector<Filament>& filaments, vec& Error,
                          double **V,  double **W, const int nt);

// let array assigment vector access pricate member vaiables
/*
   friend  void assignToFluidSolveVector(std::vector<Filament>& filaments,
                                      vec& XX, vec& YY, vec& ZZ,
                                      vec& FX, vec& FY, vec& FZ,
                                      vec& TX, vec& TY, vec& TZ);

   friend void assignToFCMarrays(const std::vector<Filament>& filaments,
                              double **Y_fcm,
                              double **F_fcm,
                              double **T_fcm);

   friend void assignToFCMarrays2(const std::vector<Filament>& filaments,
                               fcm& fcm);
 */
// same for collision barrier function that is aware of filaments
//friend void applyFilamentCollisionBarrier(std::vector<Filament>& filaments);
//friend void applyFilamentCollisionBarrierIntrinsic(std::vector<Filament>& filaments);

// same for filament property application function.
friend void apply_filament_properties(std::vector<Filament>& filaments);

// Linked list functions
friend void LinkedList::set_global_segment_number(vector<Filament>& filaments); /*,
	                                                                           int global_segment_number[Nsw][max_Nworm],
	                                                                           int global_segment_to_local_filament_number[Nsw*max_Nworm],
	                                                                           int global_segment_to_local_segment_number[Nsw*max_Nworm]);*/

friend void LinkedList::link(vector<Filament>& filaments);

// CONSTRUCTORS
// default constructor
Filament()
	: myStrainTwist(3)
	, myNworm(Nworm)
	, L_filament((Nworm - 1)*dL)
	, L_filament_long(Nworm*dL)
	, myKAP(KAP)
	, AngFrequency(KAP*(Sp4)/(4*PI*mu*pow((Nworm*dL),4.))) // omega
	, Phase(0.)
	, Wavenumber(2*K0*PI/(Nworm*dL))
	, Alpha(SwimmingHelixAlpha)
	, Beta(SwimmingHelixBeta)
	, forceLengthFilled(0)
	, stateLengthFilled(0)
	, lam(NlamPerFilament, 3,fill::zeros)
	, lam1(NlamPerFilament, 3,fill::zeros)
	, lam2(NlamPerFilament, 3,fill::zeros)
	, myNlam(NlamPerFilament)
	, myC(KAP)
	, myRhead(1.)
{
	myStrainTwist << STRAINTWISTZ << STRAINTWISTX << STRAINTWISTY;
	beads = new Bead[Nworm];
};

// constructor aware of the length
Filament(const int NumberOfBeads)
// NlamPerFilament set to NumberOfBeads-1.
	: myStrainTwist(3)
	, myNworm(NumberOfBeads)
	, L_filament((NumberOfBeads - 1)*dL)
	, L_filament_long(NumberOfBeads*dL)
	, myKAP(KAP)
	, AngFrequency(KAP*(Sp4)/(4*PI*mu*pow((NumberOfBeads*dL),4.))) // omega
	, Phase(0.)
	, Wavenumber(2*K0*PI/(NumberOfBeads*dL))
	, Alpha(SwimmingHelixAlpha)
	, Beta(SwimmingHelixBeta)
	, forceLengthFilled(0)
	, stateLengthFilled(0)
	, lam(NumberOfBeads-1, 3,fill::zeros)
	, lam1(NumberOfBeads-1, 3,fill::zeros)
	, lam2(NumberOfBeads-1, 3,fill::zeros)
	, myJapproxInv(6*NumberOfBeads, 6*NumberOfBeads,fill::zeros)
	, myNlam(NumberOfBeads-1)
	, myC(KAP)
	, myRhead(1.)
{
	myStrainTwist << STRAINTWISTZ << STRAINTWISTX << STRAINTWISTY;
	beads = new Bead[NumberOfBeads];
};

// destructor
virtual ~Filament(){
	//delete[] beads;
	// Removed because it caused segfaults when `push_back`ed. See
	// https://stackoverflow.com/questions/41712386/segmentation-fault-when-push-back-to-vector-c
}

int length() const;

Filament& print_properties() {
	cout << "myKAP" << myKAP << " ";
	cout << "myC" << myC << " ";
	cout << "Phase" << Phase << " ";
	cout << "Alpha" << Alpha << " ";
	cout << "Beta" << Beta << " ";
	cout << "forceLengthFilled" << forceLengthFilled << " ";
	cout << "stateLengthFilled" << stateLengthFilled << " ";
	cout << endl;
}

// set up all beads in an initially straight filament given head's position
// tangent vector and normal vector.
Filament& initialFilamentSetup(const vec& HeadPos,
                               const vec& tangent, const vec& normal,
                               const int forceArrayLengthFilled,
                               const int stateArrayLengthFilled);

// set up in random initial orientation given head position
Filament& initialFilamentSetupIsoRand(const vec& HeadPos,
                                      const int forceArrayLengthFilled,
                                      const int stateArrayLengthFilled,
                                      double phase_in);

// random initial orientation with the centre at the spawn point
Filament& initialFilamentSetupIsoRandCentred(const vec& HeadPos,
                                             const int forceArrayLengthFilled,
                                             const int stateArrayLengthFilled);

// random initial orientation with the centre at the spawn point
Filament& initialFilamentSetupBetweenTwoPoints(const vec& node1,
                                               const vec& node2,
                                               const int forceArrayLengthFilled,
                                               const int stateArrayLengthFilled);

Filament& copy(const Filament& other){
	// copy state variables of other filament (need to be same length!)
	for(int i=0; i<myNworm; ++i) {
		beads[i].copy(other.beads[i]);
	}
	forceLengthFilled = other.forceLengthFilled;
	stateLengthFilled = other.stateLengthFilled;
	lam = other.lam;
	return *this;
}


Filament& copyState(const Filament& other){
	// copy state variables (not Lagrange multipliers) of other filament
	assert(myNworm == other.myNworm);
	for(int i=0; i<myNworm; ++i) {
		beads[i].copy(other.beads[i]);
	}

	return *this;
}


// initial guess for time step
Filament& initialGuess(){

	// update position and orientation of first bead
	beads[0].initialGuess();

	// for all subsequent beads, update orientation only.
	for(int i = 1; i<myNworm; ++i) {
		beads[i].initialGuessRobotArm();
	}
	return *this;
}

// // initial guess for time step
// Filament& improveGuessUsingVelocity(){
//      for(int i = 0; i<myNworm; ++i) {
//              beads[i].improvedGuess();
//      }
//      return *this;
// }

// get swimmer centre of mass
vec getCOM() const {
	vec COM(3,fill::zeros);
	for(int n = 0; n<myNworm; ++n) {
		COM += beads[n].X;
	}
	COM /= myNworm;
	return COM;
}

// Check overlap with all other filaments in array otherFilaments[]
// until NumberOthers array element
bool checkOverlapWith(const std::vector<Filament>& otherFilaments,
                      const int NumberOthers,
                      float sep_factor) const;

// Check overlap with all other filaments in array otherFilaments[]
// until NumberOthers array element. Takes into account periodic boundary
// conditions. Not very efficiently written, as called only once.
bool checkOverlapWithPeriodic(const std::vector<Filament>& otherFilaments,
                              const int NumberOthers,
                              float sep_factor);

// set forces and torques acting on all beads to zero.
Filament& setZeroForcesTorques(){
	for(int i = 0; i<myNworm; ++i) {
		beads[i].F.zeros();
		beads[i].TAU.zeros();
	}
	return *this;
}

// Added for the gait reset stuff. Not normally used.
Filament& setZeroLambdas(){
	for(int i = 0; i < myNlam; ++i) {
		lam(i,0) = 0;
		lam(i,1) = 0;
		lam(i,2) = 0;
	}
}

// elastic torques from our discretisation of beam theory
Filament& applyElasticTorques(const int nt){

	UnitQuaternion q_midpoint;
	UnitQuaternion q2;
	vec vtmp(3,fill::zeros);

	double tt  = nt*dt;

	for(int i = 0; i<(myNworm-1); ++i) {
		int j = i+1;

		// Unit conversion of curvature
		double T0 = (curvature/this->L_filament_long);

		double s = dL*i + .5*dL; // arc-length to centre of segment
		double facPref; // Curvature decay factor on Alpha & Beta (swimming numbers)
		if(SpermCurvatureDecay) {
			if (s > this->L_filament/2.0) {
				facPref = T0*2.0*(this->L_filament-s)/this->L_filament;
			} else {
				facPref = T0;
			}
		} else {
			// Note lack of 2.0 here for historical reasons
			facPref = T0*(this->L_filament-s)/this->L_filament;
		}

		q2  = beads[j].qs;
		q_midpoint.assignMidPointOf(beads[i].qs,q2);

		q2 -= beads[i].qs;
		q2 *= oneOverDL;

		q2.quaproFromLeftInPlaceConj(q_midpoint); // [q2 = (q_midpoint^*)*q2]

		// swimming via time-dependent strain-twist vector
		myStrainTwist(1) =  Beta*facPref*std::cos(Wavenumber*s - AngFrequency*tt + Phase);
		myStrainTwist(2) = Alpha*facPref*std::sin(Wavenumber*s - AngFrequency*tt + Phase);

		vtmp = 2.*q2.get3Vector() - myStrainTwist;
		//cout << "[A] " << 2.*q2.get3Vector() << endl;
		//cout << "[B] " << myStrainTwist << endl;

		vtmp(0) *= myC;
		vtmp(1) *= myKAP;
		vtmp(2) *= myKAP;

		q_midpoint.RotateVector(vtmp);

		beads[i].TAU += vtmp;
		beads[j].TAU -= vtmp;
	}
	return *this;
}

// inextensibility constraint forces and torques
Filament& applyConstraintForcesTorques(){

	double fac = -.5*dL;
	vec lam_tmp(3,fill::zeros); // Lagrange mult. vector
	vec tan_tmp(3,fill::zeros); // tangent vector

	for(int i = 0; i<myNlam; ++i) {
		int j = i+1;

		// constraint forces
		lam_tmp = (lam.row(i)).t();
		beads[i].F -= lam_tmp;
		beads[j].F += lam_tmp;

		// constraint torques

		beads[i].qs.assignTangentVec(tan_tmp); // directly computes tangent
		vec tan_tmp2 = cross(tan_tmp,lam_tmp);
		beads[i].TAU += fac*tan_tmp2; //XXX

		beads[j].qs.assignTangentVec(tan_tmp); // directly computes tangent
		tan_tmp2 = cross(tan_tmp,lam_tmp);
		beads[j].TAU += fac*tan_tmp2;
	}
	return *this;
}


// External forces and torques
// NOTE: This is no longer used for swimming forces/torques.
// Instead, set SwimmingHelixAlpha and SwimmingHelixBeta to nonzero,
// and the driving will be done by applyElasticTorques (the elastic
// forces/torques will try to make the worm into a specific shape at
// a given time, which is equivalent).
Filament& applyExternalForcesTorques(const int nt){

	// for gravity:
	#if SedimentationProblem
		for(int i = 0; i<myNworm; ++i) {
			beads[i].F(0) += WeightPerLengthX*DL;
			beads[i].F(1) += WeightPerLengthY*DL;
			beads[i].F(2) += WeightPerLengthZ*DL;
		}
	#endif

	// for gravity on the swimmer only
	#if FrozenSwimmer
		if(Alpha > 1e-6 || Beta > 1e-6) {
			for(int i = 0; i<myNworm; ++i) {
				beads[i].F(0) += WeightPerLengthX*DL;
				beads[i].F(1) += WeightPerLengthY*DL;
				beads[i].F(2) += WeightPerLengthZ*DL;
			}
		}
	#endif

	return *this;
}


Filament& applyTestForcesTorques(const int nt){

	return *this;
}


void getOutOfStateVec(const vec& x, vec& xPart) const {
	// copy the part of state vector x corresponding to this filament into xPart.
	xPart = x.subvec(stateLengthFilled, stateLengthFilled+6*myNworm-1);
	return;
}


void putIntoStateVec(vec& x, const vec& xPart) const {
	// copy xPart into the part of state vector x corresponding to this filament.
	x.subvec(stateLengthFilled, stateLengthFilled+6*myNworm-1) = xPart;
	return;
}


void update(const vec& uState) {
	// update filament state by state vector uState

	vec Xupdate(3);
	vec Uupdate(3);

	// first bead gets update in position and LieAlgebra element
	int j = stateLengthFilled; // + 0 (0th index in state array)
	Xupdate << uState(j)   << uState(j+1) << uState(j+2);
	Uupdate << uState(j+3) << uState(j+4) << uState(j+5);
	beads[0].updateWith(Xupdate,Uupdate);

	// lie algebra update for remaining beads
	j = stateLengthFilled + 6;
	for(int i = 1; i < myNworm; ++i) {
		Uupdate << uState(j) << uState(j+1) << uState(j+2);
		beads[i].updateWithLie(Uupdate);
		j += 3;
	}

	for(int i = 0; i < myNlam; ++i) {
		lam(i,0) -= uState(j);
		lam(i,1) -= uState(j+1);
		lam(i,2) -= uState(j+2);
		j += 3;
	}
	return;
}


Filament& step(int nt) {
	// shift current state variables to past state variables.
	for(int i = 0; i<myNworm; ++i) {
		beads[i].step();
	}

	// shift current lagrange multipliers to past lagrange multipliers.
	if(nt>12) {
		vec lamtmp  = 3*lam - 3*lam1 + lam2;
		lam2 = lam1;
		lam1 = lam;
		lam  = lamtmp;
	}
	else if(nt==12) {
		vec lamtmp = 2.*lam - lam1;
		lam2 = lam1;
		lam1 = lam;
		lam  = lamtmp;
	}
	else {
		lam1 = lam;
	}
	return *this;
}


Filament& updateQuaternions(){
	for(int i = 0; i<myNworm; ++i) {
		beads[i].updateQuaternion();
	}
	return *this;
}


Filament& unRobotArmify(){
	// get positions of beads from head position and all orientations.
	vec Xsum(3);
	vec tangent1(3);
	vec tangent2(3);

	Xsum = beads[0].Xs;
	beads[0].qs.assignTangentVec(tangent1);

	for(int i = 1; i < myNworm; ++i) {
		beads[i].qs.assignTangentVec(tangent2);
		Xsum += 0.5*dL*(tangent1 + tangent2);
		beads[i].Xs = Xsum;
		tangent1 = tangent2;
	}
	return *this;
}


// write variable names to file (first row)
void writeDataNames(ofstream& OutputFile) const {
	for (int i=0; i<myNworm; ++i) {
		OutputFile << "X" << " " << "Y" << " " <<
		        "Z" << " " << "qr" << " " << "qi" << " " <<
		        "qj" << " " << "qk" << " ";
	}
	return;
}

// write variable names to file (first row)
void writeVelocityDataNames(ofstream& OutputFile) const {
	for (int i=0; i<myNworm; ++i) {
		OutputFile << "Vx" << " " << "Vy" << " " <<
		        "Vz" << " " << "Wx" << " " << "Wy" << " " <<
		        "Wz" << " " << "Fx" << " " << "Fy" << " " << "Fz" << " ";
	}
	return;
}

// write variable names to file (first row)
void writeCollisionDataNames(ofstream& OutputFile) const {
	for (int i=0; i<myNworm; ++i) {
		OutputFile << "Fx" << " " << "Fy" << " " << "Fz" << " ";
	}
	return;
}

// write variable names to file (first row)
void writeForceDataNames(ofstream& OutputFile) const {
	for (int i=0; i<myNworm; ++i) {
		OutputFile << "Fx" << " " << "Fy" << " " << "Fz" << " "
		           << "Tx" << " " << "Ty" << " " << "Tz" << " ";
	}
	return;
}


// write to file
void writeData(ofstream& OutputFile) const {
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	for (int i=0; i<myNworm; ++i) {
		OutputFile << beads[i].Xs(0) << " "
		           << beads[i].Xs(1) << " "
		           << beads[i].Xs(2) << " "
		           << beads[i].getQuatComponent(0) << " "
		           << beads[i].getQuatComponent(1) << " "
		           << beads[i].getQuatComponent(2) << " "
		           << beads[i].getQuatComponent(3) << " ";
	}
	// OutputFile << endl;
	return;
}


// write to file
void writeVelocityData(ofstream& OutputFile, double **V, double **W, double **F, int nt) const {
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	int filament_id = forceLengthFilled/6; // Check - not sure if this is right
	if(nt!=0) {
		for (int i=0; i<myNworm; ++i) {
			OutputFile << V[filament_id][3*i] << " "
			           << V[filament_id][3*i+1] << " "
			           << V[filament_id][3*i+2] << " "
			           << W[filament_id][3*i] << " "
			           << W[filament_id][3*i+1] << " "
			           << W[filament_id][3*i+2] << " "
			           << F[filament_id][3*i] << " "
			           << F[filament_id][3*i+1] << " "
			           << F[filament_id][3*i+2] << " ";
		}
	} else {
		for (int i=0; i<myNworm; ++i) {
			OutputFile << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " ";
		}
	}

	// OutputFile << endl;
	return;
}


// write to file
void writeCollisionData(ofstream& OutputFile, double F[][3], int nt) const {
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	int filament_id = forceLengthFilled/6; // Check - not sure if this is right
	if(nt!=0) {
		for (int i=0; i<myNworm; ++i) {
			OutputFile << F[max_Nworm*filament_id + i][0] << " "
			           << F[max_Nworm*filament_id + i][1] << " "
			           << F[max_Nworm*filament_id + i][2] << " ";
		}
	} else {
		for (int i=0; i<myNworm; ++i) {
			OutputFile << 0 << " "
			           << 0 << " "
			           << 0 << " ";
		}
	}

	// OutputFile << endl;
	return;
}

// write to file
void writeForceData(ofstream& OutputFile, double **F, double **T, int nt) const {
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	int filament_id = forceLengthFilled/6; // Check - not sure if this is right
	if(nt!=0) {
		for (int i=0; i<myNworm; ++i) {
            OutputFile << F[filament_id][3*i] << " "
			           << F[filament_id][3*i+1] << " "
			           << F[filament_id][3*i+2] << " "
			           << T[filament_id][3*i] << " "
			           << T[filament_id][3*i+1] << " "
			           << T[filament_id][3*i+2] << " ";
		}
	} else {
		for (int i=0; i<myNworm; ++i) {
			OutputFile << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " "
			           << 0 << " ";
		}
	}

	// OutputFile << endl;
	return;
}


// write X, Xt data to file. By this point, Xs=X, so we back up X and Xt.
void writeBackupData(ofstream& OutputFile) const {
	OutputFile.precision(16);
	OutputFile.setf(ios::fixed);
	OutputFile.setf(ios::showpoint);
	for (int i=0; i<myNworm; ++i) {
		OutputFile << beads[i].X(0) << " " // Current timestep
		           << beads[i].X(1) << " "
		           << beads[i].X(2) << " "
		           << beads[i].Xt(0) << " " // Previous timestep
		           << beads[i].Xt(1) << " "
		           << beads[i].Xt(2) << " "
		           << beads[i].getQuatComponent(0) << " "
		           << beads[i].getQuatComponent(1) << " "
		           << beads[i].getQuatComponent(2) << " "
		           << beads[i].getQuatComponent(3) << " "
		           << beads[i].getQuatComponentQt(0) << " "
		           << beads[i].getQuatComponentQt(1) << " "
		           << beads[i].getQuatComponentQt(2) << " "
		           << beads[i].getQuatComponentQt(3) << " "
		           << beads[i].U(0) << " "
		           << beads[i].U(1) << " "
		           << beads[i].U(2) << " "
		           << beads[i].Ut(0) << " "
		           << beads[i].Ut(1) << " "
		           << beads[i].Ut(2) << " ";
	}
	for (int i=0; i<myNlam; ++i) {
		OutputFile << lam (i,0) << " " << lam (i,1) << " " << lam (i,2) << " "
		           << lam1(i,0) << " " << lam1(i,1) << " " << lam1(i,2) << " "
		           << lam2(i,0) << " " << lam2(i,1) << " " << lam2(i,2) << " ";
	}
	return;
}


// Read in X, Xt data from file. By this point, Xs=X, so we back up X and Xt.
void continue_from_checkpoint(ifstream& input_file, int forceArrayLengthFilled,
                              int stateArrayLengthFilled) {
	string temp,temp2,temp3,temp4;
	for (int i=0; i<myNworm; ++i) {
		beads[i].radius = a;
		getline(input_file,temp,' ');
		beads[i].X(0) = stod(temp); // stod: string to double
		getline(input_file,temp,' ');
		beads[i].X(1) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].X(2) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].Xt(0) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].Xt(1) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].Xt(2) = stod(temp);
		beads[i].Xs = beads[i].X;
		getline(input_file,temp,' ');
		getline(input_file,temp2,' ');
		getline(input_file,temp3,' ');
		getline(input_file,temp4,' ');
		beads[i].qs = UnitQuaternion(stod(temp),stod(temp2),stod(temp3),stod(temp4));
		beads[i].q = beads[i].qs;
		getline(input_file,temp,' ');
		getline(input_file,temp2,' ');
		getline(input_file,temp3,' ');
		getline(input_file,temp4,' ');
		beads[i].qt = UnitQuaternion(stod(temp),stod(temp2),stod(temp3),stod(temp4));
		getline(input_file,temp,' ');
		beads[i].U(0) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].U(1) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].U(2) = stod(temp);
		beads[i].Us = beads[i].U;
		getline(input_file,temp,' ');
		beads[i].Ut(0) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].Ut(1) = stod(temp);
		getline(input_file,temp,' ');
		beads[i].Ut(2) = stod(temp);
	}
	for (int i=0; i<myNlam; ++i) {
		getline(input_file,temp,' ');
		lam (i,0) = stod(temp);
		getline(input_file,temp,' ');
		lam (i,1) = stod(temp);
		getline(input_file,temp,' ');
		lam (i,2) = stod(temp);
		getline(input_file,temp,' ');
		lam1(i,0) = stod(temp);
		getline(input_file,temp,' ');
		lam1(i,1) = stod(temp);
		getline(input_file,temp,' ');
		lam1(i,2) = stod(temp);
		getline(input_file,temp,' ');
		lam2(i,0) = stod(temp);
		getline(input_file,temp,' ');
		lam2(i,1) = stod(temp);
		getline(input_file,temp,' ');
		lam2(i,2) = stod(temp);
	}
	this->forceLengthFilled = forceArrayLengthFilled;
	this->stateLengthFilled = stateArrayLengthFilled;
}


// print for DEBUGGING
void printAllQuaternions(){
	for(int i = 0; i<myNworm; ++i) {
		beads[i].q.print();
	}
	return;
}


void printEverything(){
	for(int i = 0; i<myNworm; ++i) {
		std::cout << endl;
		std::cout << "qs ("<<i<<") ";
		beads[i].qs.printOneLine();
		std::cout << endl;
		std::cout << "Xs ("<<i<<") ";
		std::cout << "[" << beads[i].Xs(0) << ", "
		          << beads[i].Xs(1) << ", "
		          << beads[i].Xs(2) << "]" << endl;
		std::cout << "Us ("<<i<<") ";
		std::cout << "[" << beads[i].Us(0) << ", "
		          << beads[i].Us(1) << ", "
		          << beads[i].Us(2) << "]" << endl;
		std::cout << "F  ("<<i<<") ";
		std::cout << "[" << beads[i].F(0) << ", "
		          << beads[i].F(1) << ", "
		          << beads[i].F(2) << "]" << endl;
		std::cout << "TAU("<<i<<") ";
		std::cout << "[" << beads[i].TAU(0) << ", "
		          << beads[i].TAU(1) << ", "
		          << beads[i].TAU(2) << "]" << endl;
	}
	std::cout << "lambdas\n";
	std::cout << lam<< endl;
	std::cout << "L_filament\n";
	std::cout << L_filament<< endl;
	return;
}


// compute explicit Jacobian for single filament given current state.
void buildApproximateJacobian(int nt){

	int Nvars = 6*myNworm;
	mat J(Nvars,Nvars,fill::zeros);

	// compute explicit/analytical approximate Jacobian J

	// We start with the entries of the Jacobian corresponding to the update equation for the position of the first particle.

	int j = 3*myNworm + 3;
	double temp = dt/(9*PI*MU*(this->beads[0].radius));

	J(0,0) = 1;
	J(1,1) = 1;
	J(2,2) = 1;

	// J(0,j)   = temp;
	// J(1,j+1) = temp;
	// J(2,j+2) = temp;

	double Jprefac = 10.; // pre-factor for specific Jacobian entries to speed up convergence

	J(0,j)   = temp/dt;// /Jprefac;
	J(1,j+1) = temp/dt;// /Jprefac;
	J(2,j+2) = temp/dt;// /Jprefac;

	if(myNworm > 1) {
		// terms related to the derivative of the constraints with respect to the Lie algebra elements.
		J(3,3,size(3*myNworm-3,3*myNworm)) = BlockC(0.5*DL);
		// J(3,3,size(3*myNworm-3,3*myNworm)) = BlockC(-0.75*DL/dt);

		// block relating the constraints to the Lagrange multipliers.
		// J(3,3*myNworm+3,size(3*myNworm-3,3*myNworm-3)) = BlockD(nt);
		J(3,3*myNworm+3,size(3*myNworm-3,3*myNworm-3)) = Jprefac*BlockD(nt)/dt;
	}

	// block encoding the dependence of the Lie algebra update equations on the Lie algebra elements.
	if(myNworm > 1) {
		J(3*myNworm,3,size(3*myNworm,3*myNworm)) = BlockE();
	} else {
		mat E(3,3,fill::zeros);
		E(0,0) = 1;
		E(1,1) = 1;
		E(2,2) = 1;
		J(3,3,size(3,3)) = E;
	}

	// derivatives of the Lie algebra update equations with respect to the Lagrange multipliers.
	if(myNworm > 1) {
		// J(3*myNworm,3*myNworm+3,size(3*myNworm,3*myNworm-3)) = BlockF();
		J(3*myNworm,3*myNworm+3,size(3*myNworm,3*myNworm-3)) = BlockF()/dt;
	}

	// cx_vec eigval;
	// cx_mat eigvec;
	// eig_gen(eigval, eigvec, J);
	// std::cout << "J eigenvalues: " << endl;
	// std::cout << eigval << endl;

	this->myJapproxInv = inv(J);
	if(myNworm > 1) {
		this->myJapproxInv(3*myNworm + 3,0,size(3*myNworm-3,6*myNworm)) /= dt; // revert pre-conditioning
	}
	return;
}


mat rcross(const mat& xx) const { // this need not be a member function (but I'm lazy).
	// Given a 3-vector x, this function produces the matrix A such that Av = cross(v,x) for any 3-vector v.
	// N.B. The input type is given as mat so that both row and column vectors work.
	// needs C++14 stdandard library.

	mat A = {
		{0., xx(2), -xx(1)},
		{-xx(2), 0., xx(0)},
		{xx(1), -xx(0), 0.}
	};
	return A;
}


// the only part that changes is the input factor in front of the tangent vectors.
mat BlockC(const double fac) const {
	mat Cc(3*myNworm - 3, 3*myNworm, fill::zeros);

	vec tk(3);
	vec tkm1(3);

	this->beads[1].qs.assignTangentVec(tk);
	this->beads[0].qs.assignTangentVec(tkm1);

	Cc(0,0,size(3,6)) = fac*join_horiz(rcross(tkm1),rcross(tk));

	for (int k=2; k<myNworm; k++) {
		int i = 3*(k-1);
		Cc(i,0,size(3,i+3)) = Cc(i-3,0,size(3,i+3));
		this->beads[k].qs.assignTangentVec(tk);
		this->beads[k-1].qs.assignTangentVec(tkm1);
		Cc(i,i,size(3,6)) += fac*join_horiz(rcross(tkm1),rcross(tk));
	}

	// std::cout << "BlockC:" << endl;
	// std::cout << Cc << endl;
	// std::cin.get();

	return Cc;
}


mat BlockD(int nt) const {
	mat D(3*myNworm-3,3*myNworm-3,fill::zeros);

	// const double fac1 = -2.*dt/(18.*PI*MU*(this->beads[0].radius));
	// if(nt <= ImplEulerSteps) { // implicit euler steps
	// fac1 *= (3./2.);
	// }
	for (int n = 0; n<myNworm-1; ++n) {
		// double facn = 1/(6*PI*MU*(this->beads[n].radius));
		double facn = -2.*dt/(18.*PI*MU*(this->beads[n].radius));
		if(nt <= ImplEulerSteps) { // implicit euler steps
			facn *= (3./2.);
		}

		int i = 3*n;
		// D(i,0)   = fac1;
		// D(i+1,1) = fac1;
		// D(i+2,2) = fac1;

		D(i,i)     = facn;// In-place addition of this sub-block covers the first particle case.
		D(i+1,i+1) = facn;
		D(i+2,i+2) = facn;

		if (n<myNworm-2) {
			// double facnp1 = 1/(6*PI*MU*(this->beads[n+1].radius));
			facn = 2.*dt/(18.*PI*MU*(this->beads[n+1].radius));
			if(nt <= ImplEulerSteps) { // implicit euler steps
				facn *= (3./2.);
			}

			D(i,i+3)   = facn;
			D(i+1,i+4) = facn;
			D(i+2,i+5) = facn;
		}
	}

	return D;
}


mat BlockE() const {

	mat E(3*myNworm,3*myNworm,fill::zeros);

	const vec StrainTwist = this->myStrainTwist;

	// We deal with the end particles seperately.
	// n = 0
	double Tfac = 1/(8*PI*MU*pow((this->beads[0].radius),3));
	vec t(3);
	vec n(3);
	vec b(3);

	vec tp1(3);
	vec np1(3);
	vec bp1(3);

	this->beads[0].qs.assignTangentVec(t);
	this->beads[0].qs.assignNormalVec(n);
	this->beads[0].qs.assignBinormalVec(b);

	this->beads[1].qs.assignTangentVec(tp1);
	this->beads[1].qs.assignNormalVec(np1);
	this->beads[1].qs.assignBinormalVec(bp1);

	// Diagonal block
	rowvec lamvec = this->lam.row(0);
	mat ConstraintPart = -0.5*DL*rcross(lamvec)*rcross(t);

	double beta = -StrainTwist(2) + 0.5*(np1(0)*b(0) + np1(1)*b(1) + np1(2)*b(2)
	                                     - n(0)*bp1(0) - n(1)*bp1(1) - n(2)*bp1(2))/DL;

	rowvec CrossPart = {b(1)*np1(2) - b(2)*np1(1) - n(1)*bp1(2) + n(2)*bp1(1),  b(2)*np1(0) - b(0)*np1(2) - n(2)*bp1(0) + n(0)*bp1(2),  b(0)*np1(1) - b(1)*np1(0) - n(0)*bp1(1) + n(1)*bp1(0)};

	mat ElasticPart = 0.5*myC*(0.5*(t + tp1)*CrossPart/DL + beta*rcross(t));
	ElasticPart += 0.5*myKAP*(rcross(tp1)*rcross(t)/DL - 0.5*(StrainTwist(0)*rcross(n) + StrainTwist(1)*rcross(b)));
	mat dOmega = Tfac*(ElasticPart + ConstraintPart);

	mat LieMat = rcross((this->beads[0].Us)).t();
	mat OmegaMat = rcross(Tfac*(this->beads[0].TAU));
	mat Block = -2*dt*(dOmega - 0.5*(OmegaMat + LieMat*dOmega) + (rcross(OmegaMat*this->beads[0].Us) + LieMat*OmegaMat + LieMat*LieMat*dOmega)/12.)/3.;

	E(0,0,size(3,3)) = Block;

	// Off-diagonal block
	ElasticPart = 0.5*myC*(beta*rcross(tp1) - 0.5*(t + tp1)*CrossPart/DL);
	ElasticPart += 0.5*myKAP*(rcross(t).t()*rcross(tp1)/DL - 0.5*(StrainTwist(0)*rcross(np1) + StrainTwist(1)*rcross(bp1)));
	dOmega = Tfac*ElasticPart;
	Block = -2*dt*(dOmega - 0.5*LieMat*dOmega + LieMat*LieMat*dOmega/12.)/3.;
	E(0,3,size(3,3)) = Block;

	// End particle
	Tfac = 1/(8*PI*MU*pow((this->beads[myNworm-1].radius),3));

	vec tm1(3);
	vec nm1(3);
	vec bm1(3);

	this->beads[myNworm-1].qs.assignTangentVec(t);
	this->beads[myNworm-1].qs.assignNormalVec(n);
	this->beads[myNworm-1].qs.assignBinormalVec(b);

	this->beads[myNworm-2].qs.assignTangentVec(tm1);
	this->beads[myNworm-2].qs.assignNormalVec(nm1);
	this->beads[myNworm-2].qs.assignBinormalVec(bm1);

	// Diagonal block
	lamvec = this->lam.row(myNworm-2);
	ConstraintPart = -0.5*DL*rcross(lamvec)*rcross(t);
	beta = -StrainTwist(2) + 0.5*(n(0)*bm1(0) + n(1)*bm1(1) + n(2)*bm1(2)
	                              - nm1(0)*b(0) - nm1(1)*b(1) - nm1(2)*b(2))/DL;
	CrossPart = {bm1(1)*n(2) - bm1(2)*n(1) - nm1(1)*b(2) + nm1(2)*b(1), bm1(2)*n(0) - bm1(0)*n(2) - nm1(2)*b(0) + nm1(0)*b(2), bm1(0)*n(1) - bm1(1)*n(0) - nm1(0)*b(1) + nm1(1)*b(0)};

	ElasticPart = -0.5*myC*(beta*rcross(t) - 0.5*(t + tm1)*CrossPart/DL);
	ElasticPart -= 0.5*myKAP*(rcross(tm1).t()*rcross(t)/DL - 0.5*(StrainTwist(0)*rcross(n) + StrainTwist(1)*rcross(b)));
	dOmega = Tfac*(ElasticPart + ConstraintPart);
	LieMat = rcross(this->beads[myNworm-1].Us).t();
	OmegaMat = rcross(Tfac*(this->beads[myNworm-1].TAU));
	Block = -2*dt*(dOmega - 0.5*(OmegaMat + LieMat*dOmega) +
	               (rcross(OmegaMat*(this->beads[myNworm-1].Us)) + LieMat*OmegaMat + LieMat*LieMat*dOmega)/12)/3;
	E(3*myNworm-3,3*myNworm-3,size(3,3)) = Block;

	// Off-diagonal block
	ElasticPart = -0.5*myC*(beta*rcross(tm1) - 0.5*(t + tm1)*CrossPart/DL);
	ElasticPart -= 0.5*myKAP*(rcross(t)*rcross(tm1)/DL - 0.5*(StrainTwist(0)*rcross(nm1) + StrainTwist(1)*rcross(bm1)));
	dOmega = Tfac*ElasticPart;
	Block = -2.*dt*(dOmega - 0.5*LieMat*dOmega + LieMat*LieMat*dOmega/12.)/3.;

	E(3*myNworm-3,3*myNworm-6,size(3,3)) = Block;
	// Now we loop over the interior particles, which experience elastic and constraint
	//  torques on both sides and thus have the most general form to their entries.
	for (int k=1; k<myNworm-1; k++) {
		Tfac = 1/(8*PI*MU*pow((this->beads[k].radius),3.));

		this->beads[k-1].qs.assignTangentVec(tm1);
		this->beads[k-1].qs.assignNormalVec(nm1);
		this->beads[k-1].qs.assignBinormalVec(bm1);

		this->beads[k].qs.assignTangentVec(t);
		this->beads[k].qs.assignNormalVec(n);
		this->beads[k].qs.assignBinormalVec(b);

		this->beads[k+1].qs.assignTangentVec(tp1);
		this->beads[k+1].qs.assignNormalVec(np1);
		this->beads[k+1].qs.assignBinormalVec(bp1);

		LieMat = rcross(this->beads[k].Us).t();
		OmegaMat = rcross(Tfac*(this->beads[k].TAU));

		// Diagonal block
		lamvec = this->lam.row(k-1) + this->lam.row(k);
		ConstraintPart = -0.5*DL*rcross(lamvec)*rcross(t);
		beta = -StrainTwist(2) + 0.5*(np1(0)*b(0) + np1(1)*b(1) + np1(2)*b(2) - n(0)*bp1(0) - n(1)*bp1(1) - n(2)*bp1(2))/DL;
		CrossPart = {b(1)*np1(2) - b(2)*np1(1) - n(1)*bp1(2) + n(2)*bp1(1), b(2)*np1(0) - b(0)*np1(2) - n(2)*bp1(0) + n(0)*bp1(2), b(0)*np1(1) - b(1)*np1(0) - n(0)*bp1(1) + n(1)*bp1(0)};
		ElasticPart = 0.5*myC*(0.5*(t + tp1)*CrossPart/DL + beta*rcross(t));
		ElasticPart += 0.5*myKAP*(rcross(tp1)*rcross(t)/DL - 0.5*(StrainTwist(0)*rcross(n) + StrainTwist(1)*rcross(b)));
		beta = -StrainTwist(2) + 0.5*(n(0)*bm1(0) + n(1)*bm1(1) + n(2)*bm1(2) - nm1(0)*b(0) - nm1(1)*b(1) - nm1(2)*b(2))/DL;
		ElasticPart -= 0.5*myC*(beta*rcross(t) - 0.5*(t + tm1)*CrossPart/DL);
		ElasticPart -= 0.5*myKAP*(rcross(tm1).t()*rcross(t)/DL - 0.5*(StrainTwist(0)*rcross(n) + StrainTwist(1)*rcross(b)));
		dOmega = Tfac*(ConstraintPart + ElasticPart);
		Block = -2*dt*(dOmega - 0.5*(OmegaMat + LieMat*dOmega) + (rcross(OmegaMat*this->beads[k].Us) + LieMat*OmegaMat + LieMat*LieMat*dOmega)/12.)/3.;
		E(3*k,3*k,size(3,3)) = Block;

		// Next, the left off-diagonal block
		ElasticPart = -0.5*myC*(beta*rcross(tm1) - 0.5*(t + tm1)*CrossPart/DL);
		ElasticPart -= 0.5*myKAP*(rcross(t)*rcross(tm1)/DL - 0.5*(StrainTwist(0)*rcross(nm1) + StrainTwist(1)*rcross(bm1)));
		dOmega = Tfac*ElasticPart;
		Block = -2*dt*(dOmega - 0.5*LieMat*dOmega + LieMat*LieMat*dOmega/12.)/3;
		E(3*k,3*k-3,size(3,3)) = Block;

		// Finally, the right off-diagonal block
		ElasticPart = 0.5*myC*(beta*rcross(tp1) - 0.5*(t + tp1)*CrossPart/DL);
		ElasticPart += 0.5*myKAP*(rcross(t).t()*rcross(tp1)/DL - 0.5*(StrainTwist(0)*rcross(np1) + StrainTwist(1)*rcross(bp1)));
		dOmega = Tfac*ElasticPart;
		Block = -2*dt*(dOmega - 0.5*LieMat*dOmega + LieMat*LieMat*dOmega/12.)/3.;
		E(3*k,3*k+3,size(3,3)) = Block;

	}

	// All that remains is to add the identity
	for (int k=0; k<3*myNworm; k++) {
		E(k,k) += 1.;
	}
	return E;
}


mat BlockF() const {
	mat F(3*myNworm,3*myNworm-3,fill::zeros);
	vec temp(3);
	for (int i=0; i<myNworm-1; i++) {
		int j = 3*i;
		double fac = dt*DL/(24*PI*MU*pow(this->beads[i].radius,3.));

		this->beads[i].qs.assignTangentVec(temp);
		mat tmat = rcross(temp).t();
		mat umat = rcross(this->beads[i].Us);

		F(j,j,size(3,3)) = fac*(tmat - 0.5*umat*tmat + umat*umat*tmat/12.);
		fac = dt*DL/(24.*PI*MU*pow(this->beads[i+1].radius,3));

		this->beads[i+1].qs.assignTangentVec(temp);
		tmat = rcross(temp).t();
		umat = rcross(this->beads[i+1].Us);
		F(j+3,j,size(3,3)) = fac*(tmat - 0.5*umat*tmat + umat*umat*tmat/12.);
	}

	return F;
}


};



#endif
