// custom unit quaternion class representing elements in SU(2) and
// algebraic operations between them, as well as operations of elements
// of its tangent space/associated Lie Algebra (su(2) ~ R^3)on them via the exponential map.

// inclusion guards
#ifndef UNITQUATNION_IS_INCLUDED
#define UNITQUATNION_IS_INCLUDED

#include "config.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

class UnitQuaternion {
// sfs14@ic.ac.uk, started 05/02018
private:

// scalar part
double scal;

//vector part
vec vect;

public:

~UnitQuaternion(){
}

// CONSTRUCTORS:

// default constructor: identity quaternion.
UnitQuaternion()
	: scal(1.)
	, vect(3,fill::zeros) {
};

// initialisation with 4-vector
// UnitQuaternion(const Eigen::Vector4d& x) : scal(x(0)), vec(x(1),x(2),x(3)) {
// };

// initialisation with 4 scalars
UnitQuaternion(const double& x0, const double& x1, const double& x2, const double& x3) : scal(x0), vect(3) {
	vect << x1 << x2 << x3;
};

// initialisation with 3-vector
UnitQuaternion(const vec& x) : scal(0.), vect(x) {
};

// initialisation with 3-vector and scalar
UnitQuaternion(const double& y, const vec& x) : scal(y), vect(x) {
};

//initialisation with quaternion
UnitQuaternion(const UnitQuaternion& p){
	scal = p.scal;
	vect = p.vect;
}

// copy quaternion into other quaternion
UnitQuaternion& copy(const UnitQuaternion& p){
	scal = p.scal;
	vect  = p.vect;
	return *this;
}

// SIMPLE OUTPUT [FOR DEBUGGING]
void print() const {
	std::cout << endl << scal << endl << vect << endl;
	return;
}

// ONE LINE OUTPUT [FOR DEBUGGING]
void printOneLine() const {
	std::cout << "[" << scal << "; " << vect(0) << ", "
									 << vect(1) << ", "
									 << vect(2) << "]";
	return;
}


// QUATERNION ARITHMETICS

// quaternion conjugation that yields new quaternion
UnitQuaternion conj() const {
	UnitQuaternion tmp = UnitQuaternion(scal, -vect);
	return tmp;
}

UnitQuaternion c() const {
	UnitQuaternion tmp = UnitQuaternion(scal, -vect);
	return tmp;
}

// replace quaternion with its conjugate
UnitQuaternion& conjInPlace()  {
	vect *= -1.;
	return *this;
}

// quaternion product that returns new quaternion ((this*p))
UnitQuaternion quapro(const UnitQuaternion& p) const {
	return UnitQuaternion(scal*p.scal - dot(vect,p.vect), scal*p.vect + p.scal*vect + cross(vect,p.vect));
}

// quaternion product that returns new quaternion ((this*p^*))
UnitQuaternion quaproConj(const UnitQuaternion& p) const {
	return UnitQuaternion(scal*p.scal + dot(vect,p.vect), -scal*p.vect + p.scal*vect - cross(vect,p.vect));
}

// quaternion product that returns new quaternion ((p*this))
UnitQuaternion quaproFromLeft(const UnitQuaternion& p) const {
	return UnitQuaternion(scal*p.scal - dot(vect,p.vect), scal*p.vect + p.scal*vect - cross(vect,p.vect));
}

// quaternion product that returns new quaternion (((p^*)*this))
UnitQuaternion quaproFromLeftConj(const UnitQuaternion& p) const {
	return UnitQuaternion(scal*p.scal + dot(vect,p.vect), -scal*p.vect + p.scal*vect + cross(vect,p.vect));
}

// inplace quaternion product (with replacement) ((this*p))
UnitQuaternion& quaproInPlace(const UnitQuaternion& p){
	double scal_tmp = scal*p.scal - dot(vect,p.vect);
	vec tmp = cross(vect,p.vect);
	vect = scal*p.vect + p.scal*vect + tmp;
	scal = scal_tmp;
	return *this;
}

// inplace quaternion product (with replacement) ((this*(p^*)))
UnitQuaternion& quaproInPlaceConj(const UnitQuaternion& p){
	double scal_tmp = scal*p.scal + dot(vect,p.vect);
	vec tmp = cross(p.vect,vect);
	vect  = -scal*p.vect + p.scal*vect + tmp;
	scal = scal_tmp;
	return *this;
}

// inplace quaternion product (with replacement) ((p*this))
UnitQuaternion& quaproFromLeftInPlace(const UnitQuaternion& p){
	double scal_tmp = scal*p.scal - dot(vect,p.vect);
	vec tmp = cross(p.vect,vect);
	vect  = scal*p.vect + p.scal*vect + tmp;
	scal = scal_tmp;
	return *this;
}

// inplace quaternion product (with replacement) (((p^*)*this))
UnitQuaternion& quaproFromLeftInPlaceConj(const UnitQuaternion& p){
	double scal_tmp = scal*p.scal + dot(vect,p.vect);
	vec tmp = cross(vect,p.vect);
	vect  = -scal*p.vect + p.scal*vect + tmp;
	scal = scal_tmp;

	return *this;
}

// CONVERSIONS TO VECTORS
vec get4Vector() const {
	vec tmp(4);
	tmp(0) = scal;
	tmp(1) = vect(0);
	tmp(2) = vect(1);
	tmp(3) = vect(2);
	return tmp;
}

void copyTo4Vector(vec& x) const {
	x(0) = scal;
	x(1) = vect(0);
	x(2) = vect(1);
	x(3) = vect(2);
	return;
}

vec get3Vector() const {
	return vect;
}

void copyTo3Vector(vec& x) const {
	x = vect;
	return;
}

// return scalar part only.
double scalar() const {
	return scal;
}

double getVecComponent(const int i) const {
	return vect(i);
}

// return norm
double qnorm() const {
	return sqrt(scal*scal + pow(norm(vect),2));
}

// euclidean inner product
double qdot(const UnitQuaternion& p) const {
	return scal*p.scal + dot(vect,p.vect);
}

// normalise (inplace)
UnitQuaternion& qnormalise() {
	double norm_tmp = scal*scal + pow(norm(vect),2);
	norm_tmp = sqrt(norm_tmp);
	// if
	if(fabs(this->qnorm()) - 1. > 1.e-12) {
		std::cout << "normalisation now \n"; // this should never be required!!
		std::cin.get();
		scal  /= norm_tmp;
		vect  /= norm_tmp;
	}
	return *this;
}

// quaternion exponential that returns new quaternion expm(v) * q
const UnitQuaternion expm(const vec& v) const {
	double theta = norm(v)*.5;
	UnitQuaternion q_update = UnitQuaternion(std::cos(theta), std::sin(theta)*normalise(v));
	return this->quaproFromLeft(q_update);
}

// inplace quaternion exponential expm(v) * q
UnitQuaternion& expmInPlace(const vec& v){
	double theta = norm(v);

	if(theta > 1.e-12) { // need large cutoff compared to machine eps.
		UnitQuaternion q_update;
		double tmp = std::sin(theta*.5)/theta;
		q_update = UnitQuaternion(std::cos(theta*.5), tmp*v);
		this->quaproFromLeftInPlace(q_update);

		return *this;

	}else{
		return *this;
	}

}

// identity quaternion
UnitQuaternion& setToIdentity(){
	scal = 1;
	vect.zeros();
	return *this;
};



// OPERATOR OVERLOAD FOR QUATERNION ARITHMETICS

// assignment  [so we can write q = p]
UnitQuaternion& operator=(const UnitQuaternion& p){
	if (this != &p) { // do not allow self-assignment
		vect  = p.vect;
		scal = p.scal;
	}
	return *this;
}

// could be improved: define += first and derive from it.

// addition [q + p]
UnitQuaternion operator+ (const UnitQuaternion& p){
	UnitQuaternion tmp(this->scal+p.scal, this->vect+p.vect);
	return tmp;
}

// subtraction, could be improved likewise
UnitQuaternion operator- (const UnitQuaternion& p){
	UnitQuaternion tmp(this->scal-p.scal, this->vect-p.vect);
	return tmp;
}

UnitQuaternion& operator-= (const UnitQuaternion& p){
	this->scal -= p.scal;
	this->vect  -= p.vect;
	return *this;
}

// scalar multiplication [q*s]
UnitQuaternion operator* (const double& s) const {
	UnitQuaternion tmp(this->scal*s, this->vect*s);
	return tmp;
}
// [s*q]
friend UnitQuaternion operator* (const double& s, const UnitQuaternion& q) {
	return q*s;
}

UnitQuaternion& operator*= (const double& s){
	this->scal *= s;
	this->vect  *= s;
	return *this;
}

// quaternion multiplocation from left q = q*p
UnitQuaternion& operator*= (const UnitQuaternion& p){
	this->quaproInPlace(p);
	return *this;
}

// quaternion multiplication [q*p]
UnitQuaternion operator* (const UnitQuaternion& p) const {
	UnitQuaternion tmp(*this);
	tmp *= p;
	return tmp;
}

// // QUATERNIONS AS ROTATIONS q -> R_q \in SO(3)
// mat RotationMatrix() const {
//      mat M(3,3);
//      std::cout << "FUNCTION NOT WRITTEN YET" << endl;
//      std::cin.get();
//      // TBC
//      return M;
// }

// directly compute R(q)*v (inplace) [replaces quaternion_rotation_mat_mult]
void RotateVector(vec& v) const {
	double vx = v(0);
	double vy = v(1);
	double vz = v(2);

	double q11 = 2*vect(0)*vect(0);
	double q22 = 2*vect(1)*vect(1);
	double q33 = 2*vect(2)*vect(2);

	// can this be made more efficient?
	v(0) = (1 - q22 - q33)*vx + 2*(vect(0)*vect(1) - scal*vect(2))*vy + 2*(vect(0)*vect(2) + scal*vect(1))*vz;
	v(1) = 2*(vect(0)*vect(1) + scal*vect(2))*vx + (1 - q11 - q33)*vy + 2*(vect(1)*vect(2) - scal*vect(0))*vz;
	v(2) = 2*(vect(0)*vect(2) - scal*vect(1))*vx + 2*(vect(1)*vect(2) + scal*vect(0))*vy + (1 - q11 - q22)*vz;

	return;
}

// directly compute R(q)*e_x  to get tangent vector (inplace)
void assignTangentVec(vec& v) const {

	double q22 = 2*vect(1)*vect(1);
	double q33 = 2*vect(2)*vect(2);

	v(0) = (1 - q22 - q33);
	v(1) = 2*(vect(0)*vect(1) + scal*vect(2));
	v(2) = 2*(vect(0)*vect(2) - scal*vect(1));

	return;
}


void assignNormalVec(vec& v) const {

	double q11 = 2*vect(0)*vect(0);
	double q33 = 2*vect(2)*vect(2);

	v(0) = 2*(vect(0)*vect(1) - scal*vect(2));
	v(1) = (1 - q11 - q33);
	v(2) = 2*(vect(1)*vect(2) + scal*vect(0));

	return;

}

void assignBinormalVec(vec& v) const {

	double q11 = 2*vect(0)*vect(0);
	double q22 = 2*vect(1)*vect(1);

	v(0) = 2*(vect(0)*vect(2) + scal*vect(1));
	v(1) = 2*(vect(1)*vect(2) - scal*vect(0));
	v(2) = (1 - q11 - q22);

	return;

}

// quaternion square root, that is, quaternion half rotation (q^.5 s.t. q^.5*q^.5 == q)
UnitQuaternion quaSquareRoot() const {
	double newscal = sqrt(.5*(1+scal));
	vec newvec  = sqrt(.5*(1-scal)) * normalise(vect);
	UnitQuaternion tmp(newscal, newvec);
	return tmp;
}

UnitQuaternion& quaSquareRootInPlace(){
	if(scal < -.9999999999999) {
		//std::cout<<endl<<"Warning! Two consecutive beads might be rotated by ~2 pi!"<<endl<<"press enter to continue"<<endl;
		//std::cin.get();
		// scal = 1.;
		// vec *= 0.;
		scal = 0.;
		//vect << 1. <<  0. << 0.;
		vect << 0. <<  0. << 1.; // Tim says this will work :-). Changed 19/3/19
		return *this;
	}
	else if(scal > .9999999999999) {
		scal = 1.;
		vect.zeros();
		return *this;
	}
	double newscal = sqrt(.5*(1+scal));
	// vect  = sqrt(.5*(1-scal)) * vect.normalise();
	vect /= sqrt(2.*(1+scal));
	scal = newscal;
	return *this;
}

UnitQuaternion& quaSquareRootTimInPlace(){
	if(scal < -.9999999999999) {
		scal = 0.;
		vect << 0.<< 0.<< 1.;
	}
	else{
		scal = sqrt(0.5*(scal+1.));
		vect /= (2.*scal);
	}
	return *this;
}


// UnitQuaternion& assignSlerpMidPointOf(const UnitQuaternion& p1, const UnitQuaternion& p2){
//      // instpired by wikipedia article on slerp.
//      double dotprod = dot(p1.vect,p2.vect); // cos of angle between quatenions
//      UnitQuaternion p2tmp(p2);
//      if (dotprod < 0.0f) {
//              p2tmp *= -1.;
//              dotprod = -dotprod;
//      }
//
//      if(fabs(p1.scal - p2.scal) < 1.e-8) {
//              this->setToIdentity();
//              return *this;
//      }
//      if(dotprod>0.9995) {
//              scal = .5*(p1.scal + p2tmp.scal);
//              vect  = .5*(p1.vect  + p2tmp.vect );
//              this->qnormalise();
//              return *this;
//      }
//
//      double theta_0 = std::cos(dotprod); // theta_0 = angle between input vectors
//      double theta = theta_0*.5; // theta = angle between p1 and result
//      double sin_theta = std::sin(theta);
//      double sin_theta_0 = std::sin(theta_0);
//
//      double s1 = std::cos(theta) - dotprod * sin_theta / sin_theta_0; // == std::sin(theta_0 - theta) / std::sin(theta_0)
//      double s2 = sin_theta / sin_theta_0;
//
//      scal = s1*p1.scal + s2*p2tmp.scal;
//      vect  = s1*p1.vect  + s2*p2tmp.vect;
//      return *this;
// }
//
// UnitQuaternion& assignLinearMidPointOf(const UnitQuaternion& p1, const UnitQuaternion& p2){
//      // double dotprod = p1.vect.dot(p2.vect); // cos of angle between quatenions
//      // UnitQuaternion p2tmp(p2);
//      // if (dot < 0.0f) {
//      // std::cout << "0.0f" << 0.0f << endl; std::cin.get();
//      // p2tmp *= -1.;
//      // dot = -dot;
//      // }
//
//      scal = .5*p1.scal + .5*p2.scal;
//      vect  = .5*p1.vect  + .5*p2.vect;
//      this->qnormalise();
//      return *this;
// }

UnitQuaternion& assignMidPointOf(const UnitQuaternion& p1, const UnitQuaternion& p2){
	*this = (p2.quaproConj(p1));
	// this -> quaSquareRootInPlace();
	this->quaSquareRootTimInPlace();
	this->quaproInPlace(p1);

  #if verbose
	if(isnan(scal)) {
		std::cout << endl << endl;
		UnitQuaternion test = p2.quaproConj(p1);
		test.print();
		std::cout << "midpoint quaternion failed\n";
		std::cin.get();
	}
  #endif

	return *this;
}

};


#endif
