#ifndef _TOY_DM_H_
#define _TOY_DM_H_

#include "Amplitudes.h"
#include "EvoEMD/EvoEMD.h"
using namespace EvoEMD;

// * Free Parameters
DECLARE_FREE_PARAMETER(Lam, 0.4);
DECLARE_FREE_PARAMETER(MX, 100);
DECLARE_FREE_PARAMETER(MH, 125);

// * All Particles
REGISTER_PARTICLE(Boson, X, 900001, 1, RETRIEVE_PARAMETER(MX), nullptr);
REGISTER_PARTICLE(Boson, H, 900025, 2, RETRIEVE_PARAMETER(MH), nullptr);

// * Particles entering the Boltzmann Equation
REGISTER_POI(900001, true);

// * Regsiter Processes
REGISTER_PROCESS(XX_HH_Amp);

#endif  //_TOY_DM_H_
