// inclusion guards
#ifndef BEAD_IS_INCLUDED
#define BEAD_IS_INCLUDED

#include "config.hpp"
#include <armadillo>
#include "UnitQuaternion.hpp"

using namespace std;
using namespace arma;

class Filament; // forward declaration

class Bead {

private:

// quaternions are private part of beads
UnitQuaternion q;
UnitQuaternion qs;
UnitQuaternion qt;

public:

// the Filament class may access quaternions
friend class Filament;

// ErrorCheck may access quaternions (e.g., for tangent vectors)
friend bool ErrorCheck(const std::vector<Filament>& filaments, vec& Error,
                       const vec VX, const vec VY, const vec VZ,
                       const vec WX, const vec WY, const vec WZ, const int nt);

vec U;
vec Us;
vec Ut;

vec X;
vec Xs;
vec Xt;

vec F;
vec TAU;

double radius;

// default constructor
Bead()
	: q(1., 0., 0., 0.)
	, qs(1., 0., 0., 0.)
	, qt(1., 0., 0., 0.)
	, U(3,fill::zeros)
	, Us(3,fill::zeros)
	, Ut(3,fill::zeros)
	, X(3,fill::zeros)
	, Xs(3,fill::zeros)
	, Xt(3,fill::zeros)
	, F(3,fill::zeros)
	, TAU(3,fill::zeros)
{
}

Bead& copy(const Bead& other){
	// copy state variables (not forces!) of other bead
	X  = other.X;
	U  = other.U;
	q  = other.q;// needed?

	Xs = other.Xs;
	Us = other.Us;
	qs = other.qs;

	Xt = other.Xt;
	Ut = other.Ut;
	qt = other.qt; // needed?

	return *this;
}

Bead& initialBeadRandom(const vec& Pos, const double& RAD) {
	X  = Pos;
	Xt = Pos;
	Xs = Pos;


	U.zeros();
	Us.zeros();
	Ut.zeros();

	radius = RAD;

	F.zeros();
	TAU.zeros();

	vec tmp(4);
	tmp.randn();
	tmp = normalise(tmp);

	UnitQuaternion ptmp(tmp(0),tmp(1),tmp(2),tmp(3)); // initialised to unity be default

	this->q  = ptmp;
	this->qs = ptmp;
	this->qt = ptmp;
	// --------------

	return *this;
}

Bead& initialBeadSetup(const vec& Pos, const vec& tangent, const vec& normal, const double& RAD) {
	X  = Pos;
	Xt = Pos;
	Xs = Pos;


	U.zeros();
	Us.zeros();
	Ut.zeros();

	F.zeros();
	TAU.zeros();

	radius = RAD;

	/* // XXX test: set all initial quaternions to unit.
	      // [the following copies the functionality of initial_quaternion_setup].
	      // quaternion that maps e_x to tangent
	      UnitQuaternion qtmp(tangent(0),0,-tangent(2),tangent(1));
	      qtmp.quaSquareRootInPlace();

	      // apply it to e_y and solve for rotation which maps e_y to normal.
	      vec vtmp(3,fill::zeros);
	      vtmp(1) = 1;
	      qtmp.RotateVector(vtmp);

	      vec crosstmp(cross(vtmp,normal));
	      double dottmp(dot(vtmp,normal));

	      UnitQuaternion ptmp(dottmp, crosstmp);
	      ptmp.quaSquareRootInPlace();
	   this->q = ptmp*qtmp;
	 */

    /* Taken out 19/3/19
	// --------------
	UnitQuaternion ptmp; // initialised to unity be default

	this->q  = ptmp;
	this->qs = ptmp;
	this->qt = ptmp;
	// --------------
    */
    // Replaced by this 19/3/19
    UnitQuaternion qtmp(tangent(0),0,-tangent(2),tangent(1));
    qtmp.quaSquareRootInPlace();
    //std::cout<<"B"<<endl;
    // apply it to e_y and solve for rotation which maps e_y to normal.
    vec vtmp(3,fill::zeros);
    vtmp(1) = 1;
    qtmp.RotateVector(vtmp);

    vec crosstmp(cross(vtmp,normal));
    double dottmp(dot(vtmp,normal));

    UnitQuaternion ptmp(dottmp, crosstmp);
    ptmp.quaSquareRootInPlace();
    this->q = ptmp*qtmp;
    this->qs = ptmp*qtmp;
	this->qt = ptmp*qtmp;
    // end replacement 19/3/19

	return *this;
}


// set up bead given initial quaternion
Bead& initialBeadSetupQuat(const vec& Pos, const UnitQuaternion& q_in, const double& RAD) {
	X  = Pos;
	Xt = Pos;
	Xs = Pos;


	U.zeros();
	Us.zeros();
	Ut.zeros();

	F.zeros();
	TAU.zeros();

	radius = RAD;

	// --------------
	UnitQuaternion ptmp(q_in);

	this->q  = ptmp;
	this->qs = ptmp;
	this->qt = ptmp;
	// --------------

	return *this;
}

// TIME EVOLUTION
Bead& initialGuess(){
	Xs = 2.*X - Xt;
	Us = 2.*U - Ut;

	qs = q;
	qs.expmInPlace(Us);

	return *this;
}

Bead& initialGuessRobotArm(){
	Us = 2.*U - Ut;
	qs = q;
	qs.expmInPlace(Us);
	return *this;
}

// Bead& improvedGuess(){
//      Xs = X + (dt/(6*PI*a*mu))*F;
//      // Xs = X + .99*(dt/(6*PI*a*mu))*F; // for convergence test
//
//      Us = (dt/(8*PI*a*a*a*mu))*TAU;
//      // Us = U + .1*(dt/(8*PI*a*a*a*mu))*TAU; // for convergence test
//
//      // std::cout << "improvement\n";
//      // std::cout << "TAU" <<endl<< TAU <<endl;
//      // std::cin.get();
//      qs = q;
//      qs.expmInPlace(Us);
//
//      return *this;
// }


Bead& updateWith(const vec& Xupdate, const vec& Uupdate){
	// update position and Lie Algebra element
	Xs -= Xupdate;
	Us -= Uupdate;
	// update quaternion.
	qs = q;
	qs.expmInPlace(Us);

	return *this;
}

Bead& updateWithLie(const vec& Uupdate){
	// updateLie Algebra element
	Us -= Uupdate;

	// update quaternion.
	qs = q;
	qs.expmInPlace(Us);

	return *this;
}

Bead& updateQuaternion(){
	// update quaternion.
	qs = q;
	qs.expmInPlace(Us);
	return *this;
}


Bead& step(){
	// update position and Lie Algebra element
	Xt = X;
	X  = Xs;
	Ut = U;
	U  = Us;
	// update quaternion.
	qt = q;
	q.expmInPlace(U);
	return *this;
}


double getQuatComponent(const int i) const {
	if(i==0) {
		return qs.scalar();
	}
	else{
		return qs.getVecComponent(i-1);
	}
}

// Previous timestep quaternion
double getQuatComponentQt(const int i) const {
	if(i==0) {
		return qt.scalar();
	}
	else{
		return qt.getVecComponent(i-1);
	}
}

};


#endif
