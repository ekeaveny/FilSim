// inclusion guard
#ifndef HEADER_IS_INCLUDED
#define HEADER_IS_INCLUDED

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#define ARMA_NO_DEBUG
#include <armadillo>
#include <math.h>
#include <rfftw_mpi.h>
//#include <fftw3-mpi.h>
#include <iomanip>

// ALL PARAMETERS SET BY CONFIG.HPP
#include "config.hpp"
#include <assert.h>

// for profiling
#include <ctime>
#include <chrono>
#include <typeinfo>

using namespace std;
using namespace arma;

// DERIVED GLOBAL SIMULATION PARAMETERS ***************************************

// constants
extern int Np; // will be changed depending on  filament initialisation
const int NlamPerFilament = (Nworm-1);
const int Nlam = (Nsw*NlamPerFilament);
const int Nbroy = (6*Np + 3*Nlam); // deprecated

// generic Filament discretisation
const double L  = (Nworm*DL);
const double L_swimmer = (Nworm_swimmer*DL);
const double oneOverDL = (1./DL);


// fcm box and parameter setup
const double sigmadip =  a/std::pow(36.0*PI, 1.0/6.0); // see eqn (32) Journal of Computational Physics 184 (2003) 381–405
const double dx = sigmadip/1.5;
const double LfcmBox_x = NPTS_X*dx;
const double LfcmBox_y = NPTS_Y*dx;
const double LfcmBox_z = NPTS_Z*dx;

const double sigmadipsq = sigmadip*sigmadip;
const double anormdip = 1.0/std::sqrt(2.0*PI*sigmadipsq);
const double anormdip2 = 2.0 * sigmadipsq;

const double sigma = a / std::sqrt(PI);
const double sigmasq = sigma * sigma;
const double anorm = 1.0/std::sqrt(2.0*PI*sigmasq);
const double anorm2 = 2.0*sigmasq;

// This is here ONLY to define dt.
// Note that omega and f are not used in the code. Instead the frequency of
// each swimmer is set in Filament.hpp (using the same formula). This allows
// swimmers of different lengths in theory, although, obviously, dt has to
// be global.
#if SwimProblem || TwoParticleProblem
    const double omega = (KAP*(Sp4)/(4*PI*mu*pow(L_swimmer,4.)));
    const double f = omega/(2*PI);
    const double dt = 1/(StepsPerPeriod*f);     /* Simulation time-step. */
#endif

#if SedimentationProblem
const double f  = 1.;
const double omega = 2.*PI;
// const double dt = L*L/(Kap*Bnumber*StepsPerSettlingTime)*10;
const double WeightPerLength = std::sqrt(WeightPerLengthX*WeightPerLengthX + WeightPerLengthY*WeightPerLengthY + WeightPerLengthZ*WeightPerLengthZ);
const double Kap = L*L*L*WeightPerLength/Bnumber;
const double C = Kap;
const double dt = MU*L/(WeightPerLength*StepsPerSettlingTime);
#endif

// ***************************************************************************

// INCLUDE CLASSES THAT DEPEND ON GLOBAL VARIABLES HERE **********************
//#include "RPYgeneric.hpp"
#include "UnitQuaternion.hpp"
#include "Bead.hpp"
#include "Filament.hpp"
// ***************************************************************************

#endif
