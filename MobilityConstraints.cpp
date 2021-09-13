#include "multi_filament_header.hpp"

// also includes ErrorCheck function (below)


// copy filament object information into force, torque position vector.
// this loops over all particles and finds the swimmer they belong to, not vice versa.
// --> serially slower, but easier to parallelise and more robust.
/*
void assignToFluidSolveVector(std::vector<Filament>& filaments,
                              vec& XX, vec& YY, vec& ZZ,
                              vec& FX, vec& FY, vec& FZ,
                              vec& TX, vec& TY, vec& TZ){

	assert(XX.size() == Np);

	// #pragma omp parallel for
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
		// std::cout << "---------\n i: " << i << endl << "---------\n";
		// ====================================================================

		// find bead

		// positions
		XX[j] = filaments[n].beads[i].Xs(0);
		YY[j] = filaments[n].beads[i].Xs(1);
		ZZ[j] = filaments[n].beads[i].Xs(2);

		// forces
		FX[j] = filaments[n].beads[i].F(0);
		FY[j] = filaments[n].beads[i].F(1);
		FZ[j] = filaments[n].beads[i].F(2);

		// torques
		TX[j] = filaments[n].beads[i].TAU(0);
		TY[j] = filaments[n].beads[i].TAU(1);
		TZ[j] = filaments[n].beads[i].TAU(2);
	}

	return;
}


void assignToFCMarrays(const std::vector<Filament>& filaments,
                       double **Y_fcm,
                       double **F_fcm,
                       double **T_fcm){

	// copies and scales torque by LfcmBox

	// assert(XX.size() == Np);

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
		// std::cout << "---------\n i: " << i << endl << "---------\n";
		// ====================================================================

		// find bead
		// positions (in box of length 2 pi, hence modulo 2 pi)

		Y_fcm[j][0] = std::fmod(filaments[n].beads[i].Xs(0),LfcmBox_x); // /LfcmBox
		Y_fcm[j][1] = std::fmod(filaments[n].beads[i].Xs(1),LfcmBox_y);
		Y_fcm[j][2] = std::fmod(filaments[n].beads[i].Xs(2),LfcmBox_z);

		// without modulo
		// Y_fcm[j][0] = filaments[n].beads[i].Xs(0); // /LfcmBox
		// Y_fcm[j][1] = filaments[n].beads[i].Xs(1);
		// Y_fcm[j][2] = filaments[n].beads[i].Xs(2);

		// forces
		F_fcm[j][0] = filaments[n].beads[i].F(0);
		F_fcm[j][1] = filaments[n].beads[i].F(1);
		F_fcm[j][2] = filaments[n].beads[i].F(2);

		// torques
		T_fcm[j][0] = filaments[n].beads[i].TAU(0);
		T_fcm[j][1] = filaments[n].beads[i].TAU(1);
		T_fcm[j][2] = filaments[n].beads[i].TAU(2);
	}

	return;
}
*/

// not part of the unit quaternion class, but related to quaternions:
vec dexpinv(const vec& u, const vec& w){
	double tmp = norm(u); //  = fmod(u.norm(),PI);
	vec out(3,fill::zeros);

	if (tmp >= 1.e-12) {
		double tmp2;
		tmp2  = tmp/std::sin(.5*tmp);// for small values, close to 1. OK
		tmp2 *= std::cos(.5*tmp); // for small values close to 0. OK
		tmp   = .5*tmp2 - 1.;// usually not subtracting similar numbers

		out = w - .5*cross(u,w);
		out -= tmp*(cross(normalise(u),cross(normalise(u),w) ));
    #if verbose
		std::cout << "dexpinv\n"; // XXX
		std::cout << out << endl; // XXX
    #endif
		return out;
	}
	else {
		// #if verbose
		// std::cout << "small u dexpinv approximation \n"; std::cin.get();
		// #endif
		// tmp *= tmp;
		tmp = 1./12;

		out = w - .5*cross(u,w);
		out  += tmp*(cross(u,cross(u,w)));
    #if verbose
		std::cout << "dexpinv (taylor)\n"; // XXX
		std::cout << out << endl; // XXX
    #endif
		return out;

	}
}





// this function contains implicit time integration with constraints.
/*
bool ErrorCheck(const std::vector<Filament>& filaments, vec& Error,
                const vec VX, const vec VY, const vec VZ,
                const vec WX, const vec WY, const vec WZ, const int nt)
{
    // Only called from rank 0.

	vec V(3);
	// vec Xacc(3,fill::zeros);
	vec W(3);
	vec U(3);
	vec X(3);
	vec Us(3);
	vec Xs(3);
	vec Err1(3);
	vec Err2(3);
	// vec tangent1(3,fill::zeros);
	// vec tangent2(3,fill::zeros);
	// vec tangent3(3,fill::zeros);
	// vec tangentBDF(3,fill::zeros);
	UnitQuaternion q;
	Bead bead_i;
	bool check = false;

	// check time step
	if(nt <= ImplEulerSteps) { // use implicit Euler for initial steps

		// const double fourthirds  = 4./3.;
		// const double onethird    = 1./3.;
		// const double twothirdsDt = dt*2./3.;
		// #pragma omp parallel for private(V,W,U,X,Us,Xs,Err1,Err2,tangent1,tangent2,q,bead_i)
		for(int nn = 0; nn<filaments.size(); ++nn) {

			// filament_n.copyState(filaments[nn]);
			int inVfirst = filaments[nn].forceLengthFilled;
			int inSfirst = filaments[nn].stateLengthFilled;
			int Nw       = filaments[nn].length();

			int j = inVfirst; // index for force/velocity vectors in [0,Nw]
			int k = inSfirst; // index for state vectors in [0,6*Nw]
			int i = 0; // index for beads in [0,Nw]

			V << VX(j) << VY(j) << VZ(j);
			W << WX(j) << WY(j) << WZ(j);

			// Xacc = bead_i.Xs;

			bead_i.copy(filaments[nn].beads[i]);

			// position impl. euler
			Err1 = bead_i.Xs - bead_i.X - dt*V;

			// Lie Algebra impl. euler
			// Err2 = bead_i.Us  - dt*dexpinv(bead_i.Us, W);
			Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);


			// // position BDF2
			// Err1 = bead_i.Xs - fourthirds*bead_i.X + onethird*bead_i.Xt - twothirdsDt*V;
			//
			// // Lie Algebra BDF2
			// // Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);
			// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);

			// first bead position update error
			Error(k+0) = Err1(0);
			Error(k+1) = Err1(1);
			Error(k+2) = Err1(2);

			// firsts bead Lie algebra update error
			Error(k+0+ 3*Nw) = Err2(0);
			Error(k+1+ 3*Nw) = Err2(1);
			Error(k+2+ 3*Nw) = Err2(2);


			bool checktmp = any(abs(Err1) > TOL);
			checktmp = (checktmp || any(abs(Err2) > TOL));

			if (checktmp) {
				// #pragma omp atomic write // this ensures that two threads do not try to write at the same time.
				check = true; // shared. if any thread writes, value changes. right? XXX
			}

			// bead_i.qs.assignTangentVec(tangent1);

			// loop over remaining beads and constraints (assuming NlamPerFilament = Nworm - 1)
			for(i = 1; i<Nw; ++i) {

				j = inVfirst+i;
				k = inSfirst+3*i;

				V << VX(j) << VY(j) << VZ(j);
				W << WX(j) << WY(j) << WZ(j);

				// tangentBDF.zeros();

				// past bead (can be optimised away!!)
				// bead_i.qs.assignTangentVec(tangent1); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula

				// Xacc += .5*dL*tangent1;


				// this bead
				bead_i.copy(filaments[nn].beads[i]);

				// bead_i.qs.assignTangentVec(tangent2); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula

				// Xacc += .5*dL*tangent2;


				// Vacc += (.75*dL/dt)*tangentBDF;

				// position impl. euler
				// Err1 = bead_i.Xs - bead_i.X - dt*V;


				// position constraint error
				Err1 = bead_i.Xs - bead_i.X - dt*V;//V - Vacc;x

				// Lie Algebra impl. euler
				Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);

				// // Lie Algebra BDF2
				// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);

				// velocity constraint error
				Error(k+0) = Err1(0);
				Error(k+1) = Err1(1);
				Error(k+2) = Err1(2);

				//Lie algebra update error
				Error(k+0+ 3*Nw) = Err2(0);
				Error(k+1+ 3*Nw) = Err2(1);
				Error(k+2+ 3*Nw) = Err2(2);


				checktmp = any(abs(Err1) > TOL);
				checktmp = (checktmp || any(abs(Err2) > TOL));

				// checktmp = (checktmp || any(abs(Err1) > TOL2));

				if (checktmp) {
					// #pragma omp atomic write
					check = true;
				}

				// tangent1 = tangent2; // avoids computing tangent vectors twice.

			}
		}
	}
	else{

		const double fourthirds  = 4./3.;
		const double onethird    = 1./3.;
		const double twothirdsDt = dt*2./3.;

		// #pragma omp parallel for private(V,W,U,X,Us,Xs,Err1,Err2,tangent1,tangent2,q,bead_i)
		for(int nn = 0; nn<filaments.size(); ++nn) {

			// filament_n.copyState(filaments[nn]);
			int inVfirst = filaments[nn].forceLengthFilled;
			int inSfirst = filaments[nn].stateLengthFilled;
			int Nw       = filaments[nn].length();

			int j = inVfirst; // index for force/velocity vectors in [0,Nw]
			int k = inSfirst; // index for state vectors in [0,6*Nw]
			int i = 0; // index for beads in [0,Nw]

			V << VX(j) << VY(j) << VZ(j);
			W << WX(j) << WY(j) << WZ(j);
			// Vacc = V;

			bead_i.copy(filaments[nn].beads[i]);

			// position impl. euler
			// Err1 = bead_i.Xs - bead_i.X - dt*V;

			// Lie Algebra impl. euler
			// Err2 = bead_i.Us  - dt*dexpinv(bead_i.Us, W);
			// Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);


			// position BDF2
			Err1 = bead_i.Xs - fourthirds*bead_i.X + onethird*bead_i.Xt - twothirdsDt*V;

			// Lie Algebra BDF2
			// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);
			Err2 = bead_i.Us - onethird*bead_i.U - twothirdsDt*dexpinv(bead_i.Us, W);


			// first bead position update error
			Error(k+0) = Err1(0);
			Error(k+1) = Err1(1);
			Error(k+2) = Err1(2);

			// firsts bead Lie algebra update error
			Error(k+0+ 3*Nw) = Err2(0);
			Error(k+1+ 3*Nw) = Err2(1);
			Error(k+2+ 3*Nw) = Err2(2);


			bool checktmp = any(abs(Err1) > TOL);
			checktmp = (checktmp || any(abs(Err2) > TOL));

			if (checktmp) {
				// #pragma omp atomic write // this ensures that two threads do not try to write at the same time.
				check = true; // shared. if any thread writes, value changes. right? XXX
			}

			// bead_i.qs.assignTangentVec(tangent1);

			// loop over remaining beads and constraints (assuming NlamPerFilament = Nworm - 1)
			for(i = 1; i<Nw; ++i) {

				j = inVfirst+i;
				k = inSfirst+3*i;

				V << VX(j) << VY(j) << VZ(j);
				W << WX(j) << WY(j) << WZ(j);

				// tangentBDF.zeros();

				// past bead (can be optimised away!!)
				// bead_i.qs.assignTangentVec(tangent1); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula


				// this bead
				bead_i.copy(filaments[nn].beads[i]);

				// bead_i.qs.assignTangentVec(tangent1); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula


				// Vacc += (.75*dL/dt)*tangentBDF;

				// position impl. euler
				// Err1 = bead_i.Xs - bead_i.X - dt*V;
				Err1 = bead_i.Xs - fourthirds*bead_i.X + onethird*bead_i.Xt - twothirdsDt*V;


				// velocity constraint error
				// Err1 = V - Vacc;

				// Lie Algebra impl. euler
				// Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);

				// // Lie Algebra BDF2
				// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);
				Err2 = bead_i.Us - onethird*bead_i.U - twothirdsDt*dexpinv(bead_i.Us, W);

				// velocity constraint error
				Error(k+0) = Err1(0);
				Error(k+1) = Err1(1);
				Error(k+2) = Err1(2);

				//Lie algebra update error
				Error(k+0+ 3*Nw) = Err2(0);
				Error(k+1+ 3*Nw) = Err2(1);
				Error(k+2+ 3*Nw) = Err2(2);


				checktmp = any(abs(Err1) > TOL);
				checktmp = (checktmp || any(abs(Err2) > TOL));

				// checktmp = (checktmp || any(abs(Err1) > TOL2));

				if (checktmp) {
					// #pragma omp atomic write
					check = true;
				}

				// tangent1 = tangent2; // avoids computing tangent vectors twice.

			}
		}
	}

	// std::cout << "Error.norm()" << norm(Error) << endl;

	return check;
}
*/



bool ErrorCheckFCM(const std::vector<Filament>& filaments, vec& Error,
                   double **V_fcm,  double **W_fcm, const int nt){

    // Only called from rank 0.

	// multiplies velocities by  LfcmBox and error checks

	vec V(3);
	// vec Xacc(3,fill::zeros);
	vec W(3);
	vec U(3);
	vec X(3);
	vec Us(3);
	vec Xs(3);
	vec Err1(3);
	vec Err2(3);
	// vec tangent1(3,fill::zeros);
	// vec tangent2(3,fill::zeros);
	// vec tangent3(3,fill::zeros);
	// vec tangentBDF(3,fill::zeros);
	UnitQuaternion q;
	Bead bead_i;
	bool check = false;
    bool checkexplode = false;

	// check time step
	if(nt <= ImplEulerSteps) { // use implicit Euler for initial steps

		// const double fourthirds  = 4./3.;
		// const double onethird    = 1./3.;
		// const double twothirdsDt = dt*2./3.;
		// #pragma omp parallel for private(V,W,U,X,Us,Xs,Err1,Err2,tangent1,tangent2,q,bead_i)
		for(int nn = 0; nn<filaments.size(); ++nn) {

			// filament_n.copyState(filaments[nn]);
			int inVfirst = filaments[nn].forceLengthFilled;
			int inSfirst = filaments[nn].stateLengthFilled;
			int Nw       = filaments[nn].length();

			int j = inVfirst; // index for force/velocity vectors in [0,Nw]
			int k = inSfirst; // index for state vectors in [0,6*Nw]
			int i = 0; // index for beads in [0,Nw]

			V << V_fcm[j][0] << V_fcm[j][1] << V_fcm[j][2];
			W << W_fcm[j][0] << W_fcm[j][1] << W_fcm[j][2];

			// Xacc = bead_i.Xs;

			bead_i.copy(filaments[nn].beads[i]);

			// position impl. euler
			Err1 = bead_i.Xs - bead_i.X - dt*V;

			// Lie Algebra impl. euler
			// Err2 = bead_i.Us  - dt*dexpinv(bead_i.Us, W);
			Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);

			// // position BDF2
			// Err1 = bead_i.Xs - fourthirds*bead_i.X + onethird*bead_i.Xt - twothirdsDt*V;
			//
			// // Lie Algebra BDF2
			// // Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);
			// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);

			// first bead position update error
			Error(k+0) = Err1(0);
			Error(k+1) = Err1(1);
			Error(k+2) = Err1(2);

			// firsts bead Lie algebra update error
			Error(k+0+ 3*Nw) = Err2(0);
			Error(k+1+ 3*Nw) = Err2(1);
			Error(k+2+ 3*Nw) = Err2(2);


			bool checktmp = any(abs(Err1) > TOL);
			checktmp = (checktmp || any(abs(Err2) > TOL));

            bool checkexplodetmp = any(abs(Err1) > 1e6);
            checkexplodetmp = (checkexplodetmp || any(abs(Err2) > 1e6));

			if (checktmp) {
				// #pragma omp atomic write // this ensures that two threads do not try to write at the same time.
				check = true; // shared. if any thread writes, value changes. right? XXX
			}
            if (checkexplode) {
                checkexplode = true;
            }

			// bead_i.qs.assignTangentVec(tangent1);

			// loop over remaining beads and constraints (assuming NlamPerFilament = Nworm - 1)
			for(i = 1; i<Nw; ++i) {

				j = inVfirst+i;
				k = inSfirst+3*i;

				V << V_fcm[j][0] << V_fcm[j][1] << V_fcm[j][2];
				W << W_fcm[j][0] << W_fcm[j][1] << W_fcm[j][2];

				// tangentBDF.zeros();

				// past bead (can be optimised away!!)
				// bead_i.qs.assignTangentVec(tangent1); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula

				// Xacc += .5*dL*tangent1;


				// this bead
				bead_i.copy(filaments[nn].beads[i]);

				// bead_i.qs.assignTangentVec(tangent2); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula

				// Xacc += .5*dL*tangent2;


				// Vacc += (.75*dL/dt)*tangentBDF;

				// position impl. euler
				// Err1 = bead_i.Xs - bead_i.X - dt*V;


				// position constraint error
				Err1 = bead_i.Xs - bead_i.X - dt*V;//V - Vacc;x

				// Lie Algebra impl. euler
				Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);

				// // Lie Algebra BDF2
				// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);

				// velocity constraint error
				Error(k+0) = Err1(0);
				Error(k+1) = Err1(1);
				Error(k+2) = Err1(2);

				//Lie algebra update error
				Error(k+0+ 3*Nw) = Err2(0);
				Error(k+1+ 3*Nw) = Err2(1);
				Error(k+2+ 3*Nw) = Err2(2);


				checktmp = any(abs(Err1) > TOL);
				checktmp = (checktmp || any(abs(Err2) > TOL));

				// checktmp = (checktmp || any(abs(Err1) > TOL2));

				if (checktmp) {
					// #pragma omp atomic write
					check = true;
				}

                checkexplodetmp = any(abs(Err1) > 1e6);
                checkexplodetmp = (checkexplodetmp || any(abs(Err2) > 1e6));
                if (checkexplodetmp) {
                    checkexplode = true;
                }

				// tangent1 = tangent2; // avoids computing tangent vectors twice.

			}
		}
	}
	else{

		const double fourthirds  = 4./3.;
		const double onethird    = 1./3.;
		const double twothirdsDt = dt*2./3.;

		// #pragma omp parallel for private(V,W,U,X,Us,Xs,Err1,Err2,tangent1,tangent2,q,bead_i)
		for(int nn = 0; nn<filaments.size(); ++nn) {

			// filament_n.copyState(filaments[nn]);
			int inVfirst = filaments[nn].forceLengthFilled;
			int inSfirst = filaments[nn].stateLengthFilled;
			int Nw       = filaments[nn].length();

			int j = inVfirst; // index for force/velocity vectors in [0,Nw]
			int k = inSfirst; // index for state vectors in [0,6*Nw]
			int i = 0; // index for beads in [0,Nw]

			V << V_fcm[j][0] << V_fcm[j][1] << V_fcm[j][2];
			W << W_fcm[j][0] << W_fcm[j][1] << W_fcm[j][2];
			// Vacc = V;

			bead_i.copy(filaments[nn].beads[i]);

			// position impl. euler
			// Err1 = bead_i.Xs - bead_i.X - dt*V;

			// Lie Algebra impl. euler
			// Err2 = bead_i.Us  - dt*dexpinv(bead_i.Us, W);
			// Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);


			// position BDF2
			Err1 = bead_i.Xs - fourthirds*bead_i.X + onethird*bead_i.Xt - twothirdsDt*V;

			// Lie Algebra BDF2
			// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);
			Err2 = bead_i.Us - onethird*bead_i.U - twothirdsDt*dexpinv(bead_i.Us, W);


			// first bead position update error
			Error(k+0) = Err1(0);
			Error(k+1) = Err1(1);
			Error(k+2) = Err1(2);

			// firsts bead Lie algebra update error
			Error(k+0+ 3*Nw) = Err2(0);
			Error(k+1+ 3*Nw) = Err2(1);
			Error(k+2+ 3*Nw) = Err2(2);


            bool checktmp = any(abs(Err1) > TOL);
			checktmp = (checktmp || any(abs(Err2) > TOL));

            bool checkexplodetmp = any(abs(Err1) > 1e6);
            checkexplodetmp = (checkexplodetmp || any(abs(Err2) > 1e6));

			if (checktmp) {
				// #pragma omp atomic write // this ensures that two threads do not try to write at the same time.
				check = true; // shared. if any thread writes, value changes. right? XXX
			}
            if (checkexplodetmp) {
                checkexplode = true;
            }



			// bead_i.qs.assignTangentVec(tangent1);

			// loop over remaining beads and constraints (assuming NlamPerFilament = Nworm - 1)
			for(i = 1; i<Nw; ++i) {

				j = inVfirst+i;
				k = inSfirst+3*i;

				V << V_fcm[j][0] << V_fcm[j][1] << V_fcm[j][2];
				W << W_fcm[j][0] << W_fcm[j][1] << W_fcm[j][2];

				// tangentBDF.zeros();

				// past bead (can be optimised away!!)
				// bead_i.qs.assignTangentVec(tangent1); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula


				// this bead
				bead_i.copy(filaments[nn].beads[i]);

				// bead_i.qs.assignTangentVec(tangent1); // new tangent vector
				// bead_i.q.assignTangentVec(tangent2); // old tangent vector
				// bead_i.qt.assignTangentVec(tangent3); // older tangent vector
				// tangentBDF += tangent1 - fourthirds*tangent2 + onethird*tangent3; // BDF2 formula


				// Vacc += (.75*dL/dt)*tangentBDF;

				// position impl. euler
				// Err1 = bead_i.Xs - bead_i.X - dt*V;
				Err1 = bead_i.Xs - fourthirds*bead_i.X + onethird*bead_i.Xt - twothirdsDt*V;


				// velocity constraint error
				// Err1 = V - Vacc;

				// Lie Algebra impl. euler
				// Err2 = bead_i.Us - dt*dexpinv(bead_i.Us, W);

				// // Lie Algebra BDF2
				// Err2 = bead_i.Us - onethird*bead_i.Ut - twothirdsDt*dexpinv(bead_i.Us, W);
				Err2 = bead_i.Us - onethird*bead_i.U - twothirdsDt*dexpinv(bead_i.Us, W);

				// velocity constraint error
				Error(k+0) = Err1(0);
				Error(k+1) = Err1(1);
				Error(k+2) = Err1(2);

				//Lie algebra update error
				Error(k+0+ 3*Nw) = Err2(0);
				Error(k+1+ 3*Nw) = Err2(1);
				Error(k+2+ 3*Nw) = Err2(2);


				checktmp = any(abs(Err1) > TOL);
				checktmp = (checktmp || any(abs(Err2) > TOL));

				// checktmp = (checktmp || any(abs(Err1) > TOL2));

				if (checktmp) {
					// #pragma omp atomic write
					check = true;
				}

                bool checkexplodetmp = any(abs(Err1) > 1e6);
                checkexplodetmp = (checkexplodetmp || any(abs(Err2) > 1e6));
                if (checkexplodetmp) {
                    checkexplode = true;
                }

				// tangent1 = tangent2; // avoids computing tangent vectors twice.

			}
		}
	}
    #if flush_on
        std::cout << std::setprecision(4) << "[Error.norm() = " << norm(Error) << "]";
    #endif
    if (checkexplode)
    {
        std::cout << "Prepare for explosion: max(|Error|) > 1e6" << endl;
    }

	return check;
}





// void DiagonalMobilityError(Filament filaments[], vec& v,
//                            vec& XX,  vec& YY,  vec& ZZ,
//                            vec& FX,  vec& FY,  vec& FZ,
//                            vec& TX,  vec& TY,  vec& TZ,
//                            vec& VX,  vec& VY,  vec& VZ,
//                            vec& WX,  vec& WY,  vec& WZ,
//                            const int nt, const bool ComputeNewForces)
// {
//      // v is the output error and if CumputeNewForces the input state.
//      // std::cout << "\n\n === DiagonalMobilityError === \n\n";
//      assert(XX.size()==Np);
//   #if verbose
//      for (int i=0; i<Nbroy; i++) {cout<<"v["<<i<<"] = "<<v(i)<<"\n"; };
//      std::cin.get();
//   #endif
//
//      if(ComputeNewForces) {
//              // this is the case called by gmres. it assumes that v is also the
//              // initial new state vector.
//              // std::cout << endl <<  v.size() << endl;
//
//       #pragma omp parallel for
//              for(int nn=0; nn<Nsw; ++nn) {
//                      filaments[nn].assignFromStateVector(v); // this copies state information from input "v"
//                      filaments[nn].updateQuaternions();
//                      filaments[nn].setZeroForcesTorques();
//
//                      filaments[nn].applyElasticTorques();
//                      // #if verbose
//                      // std::cout << "1 diag Elastic Torques\n";
//                      // filaments[nn].printEverything(); // std::cin.get();
//
//                      // #endif
//
//                      filaments[nn].applyConstraintForcesTorques();
//                      // #if verbose
//                      // std::cout << "1 diag Constraint Forces Torques\n";
//
//                      // filaments[nn].printEverything(); // std::cin.get();
//                      // #endif
//                      // filaments[nn].printEverything(); std::cin.get();
//                      filaments[nn].applyDrivingTorques(nt);
//                      // #if verbose
//                      // std::cout << "1 diag Driving Torques\n";
//                      // for (int i=0; i<Np; i++){cout<<"FX["<<i<<"] = "<<FX[i]<<"\n";};
//                      // for (int i=0; i<Np; i++){cout<<"FY["<<i<<"] = "<<FY[i]<<"\n";};
//                      // for (int i=0; i<Np; i++){cout<<"FZ["<<i<<"] = "<<FZ[i]<<"\n";};
//                      // for (int i=0; i<Np; i++){cout<<"TX["<<i<<"] = "<<TX[i]<<"\n";};
//                      // for (int i=0; i<Np; i++){cout<<"TY["<<i<<"] = "<<TY[i]<<"\n";};
//                      // for (int i=0; i<Np; i++){cout<<"TZ["<<i<<"] = "<<TZ[i]<<"\n";};
//
//                      // filaments[nn].printEverything(); std::cin.get();
//                      // #endif
//                      // cin.get();
//                      // filaments[nn].applyTestForcesTorques(nt);
//              }
//
//              // for(int nn=0; nn<Nsw; ++nn){
//              //   filaments[nn].assignToFluidSolveVector(XX, YY, ZZ, FX, FY, FZ, TX, TY, TZ);
//              // }
//              assignToFluidSolveVector(filaments, XX, YY, ZZ, FX, FY, FZ, TX, TY, TZ);
//
//
//              applyCollisionBarrier<vec>(XX, YY, ZZ, FX, FY, FZ);
//
//              // this is needed for new tangent vector in error check. XXX?
//              // #pragma omp parallel for
//              // for(int nn=0; nn<Nsw; ++nn){
//              // filaments[nn].updateQuaternions();
//              // }
//      }
//
//      // double c1;// = 1./(6*PI*a*mu);
//      // double c2 = 1./8/PI/(mu*a*a*a);
//      // double c2;
//
//      // VX = c1*FX;
//      // VY = c1*FY;
//      // VZ = c1*FZ;
//
//      VX = FX/(6.*PI*a*mu);
//      VY = FY/(6.*PI*a*mu);
//      VZ = FZ/(6.*PI*a*mu);
//
//      // WX = c2*TX;
//      // WY = c2*TY;
//      // WZ = c2*TZ;
//
//      WX = TX/(8.*PI*a*a*a*mu);
//      WY = TY/(8.*PI*a*a*a*mu);
//      WZ = TZ/(8.*PI*a*a*a*mu);
//
//      // #if notAllTorques
//      //   WX *= 0.;
//      //   WY *= 0.;
//      //   // WZ *= 0.;
//      // #endif
//     #if verbose
//      std::cout << "diagonal mobility velocities \n";
//      for (int i=0; i<Nworm; i++) {cout<<"FX["<<i<<"] = "<<FX(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"FY["<<i<<"] = "<<FY(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"FZ["<<i<<"] = "<<FZ(i)<<"\n"; };
//
//      for (int i=0; i<Nworm; i++) {cout<<"TX["<<i<<"] = "<<TX(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"TY["<<i<<"] = "<<TY(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"TZ["<<i<<"] = "<<TZ(i)<<"\n"; };
//
//      for (int i=0; i<Nworm; i++) {cout<<"VX["<<i<<"] = "<<VX(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"VY["<<i<<"] = "<<VY(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"VZ["<<i<<"] = "<<VZ(i)<<"\n"; };
//
//      for (int i=0; i<Nworm; i++) {cout<<"WX["<<i<<"] = "<<WX(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"WY["<<i<<"] = "<<WY(i)<<"\n"; };
//      for (int i=0; i<Nworm; i++) {cout<<"WZ["<<i<<"] = "<<WZ(i)<<"\n"; };
//      std::cin.get();
//     #endif
//
//      //bool check_tmp;
//      ErrorCheck(filaments, v, VX, VY, VZ, WX, WY, WZ, nt);
//     #if verbose
//      for (int i=0; i<Nbroy; i++) {cout<<"after errorcheck v["<<i<<"] = "<<v(i)<<"\n"; };
//      std::cin.get();
//     #endif
//
//      return;
// }
