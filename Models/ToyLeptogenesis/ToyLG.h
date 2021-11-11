#ifndef _TOY_LEPTOGENESIS_H_
#define _TOY_LEPTOGENESIS_H_

#include "Amplitudes.h"
#include "EvoEMD/EvoEMD.h"
#include "Parameters.h"
using namespace EvoEMD;

// * Free Parameters
DECLARE_FREE_PARAMETER(MN1, 1e13);
DECLARE_FREE_PARAMETER(Lam, 4e-3);
DECLARE_FREE_PARAMETER(Eps, 1e-6);

// * Parameters that depend on free parameter
REGISTER_PARAMETER(param_GammaN1);

// * Particles
REGISTER_PARTICLE(Fermion, N1, 900001, 2, RETRIEVE_PARAMETER(MN1), RETRIEVE_PARAMETER(GammaN1));
REGISTER_PARTICLE(Fermion, dL, 900011, 2 * 2, nullptr, nullptr, true);
REGISTER_PARTICLE(Boson, Phi, 25, 2, nullptr, nullptr);

// * Particles entering the Boltzmann Equation
REGISTER_POI(900001, 1);
REGISTER_POI(900011, 0);

// * Register Processes
REGISTER_PROCESS(N_LPhi_Amp);
REGISTER_PROCESS(delta_N_LPhi_Amp);

#endif  //_TOY_LEPTOGENESIS_H_
