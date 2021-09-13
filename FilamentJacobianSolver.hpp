#ifndef JACOBIAN_SOLVER_IS_INCLUDED
#define JACOBIAN_SOLVER_IS_INCLUDED

#include "config.hpp"
#include <armadillo>
#include "multi_filament_header.hpp"

// Jacobian solver that uses filament-level sub-jacobians
// started 06/07/17 by simon sfs14@ic.ac.uk

using namespace std;
using namespace arma;

class filamentJsolver {
private:

mat Cmat;
mat Dmat;
// mat Japprox; // block is now stored as part of each filament.

public:

// default constructor
filamentJsolver()
	: Cmat(6*Np,broyden_maxiter,fill::zeros)
	, Dmat(6*Np,broyden_maxiter,fill::zeros)
{
};

filamentJsolver(int myNp)
	: Cmat(6*myNp,broyden_maxiter,fill::zeros)
	, Dmat(6*myNp,broyden_maxiter,fill::zeros)
{
};

void buildJacobian(std::vector<Filament>& filaments, int nt) {
	// Only called by rank 0.

	// reset rank one updates
	Dmat.zeros();
	Cmat.zeros();
	// for each filament, build an invert the analytical approx Jacobian
	for(int nn=0; nn < filaments.size(); ++nn) {
		filaments[nn].buildApproximateJacobian(nt);
	}
	return;
}

void solve(std::vector<Filament>& filaments,
           const vec& b, vec& x, int curr_iter) const {
	// Only called by rank 0.

	// solves Jx = b for x, where J^-1 blocks are members of filamens and rank one
	// updates to J^-1 are stored in solver object's Dmat and Cmat.
	// (Bad Broyden's)

	// block-wise multiplication with inverse jacobian.
	for(int nn=0; nn < filaments.size(); ++nn) {
		vec xFil;
		vec bFil;
		filaments[nn].getOutOfStateVec(b,bFil); // get short bFil out of long b (6*Np)
		xFil = filaments[nn].myJapproxInv*bFil;
		filaments[nn].putIntoStateVec(x,xFil); // put short xFil into long x   (6*Np)
	}

	// bad Broyden's updates
	for(int i = 0; i < curr_iter; ++i) {
		vec dtmp = Dmat.col(i);
		vec ctmp = Cmat.col(i);
		x += ctmp * dot(dtmp,b);
	}
	return;
}


void addDmatCol(const vec& v, int ind){
	Dmat.col(ind) = v;
	return;
}

void addCmatCol(const vec& v, int ind){
	Cmat.col(ind) = v;
	return;
}


}; // end of class

#endif
