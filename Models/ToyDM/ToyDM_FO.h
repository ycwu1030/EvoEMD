#ifndef _TOY_DM_H_s
#define _TOY_DM_H_s

#include "Amplitudes.h"
#include "EvoEMD/EvoEMD.h"
using namespace EvoEMD;

// * Free Parameters
DECLARE_FREE_PARAMETER(Lam, 0.8);
DECLARE_FREE_PARAMETER(MX, 100);

// * All Particles
REGISTER_PARTICLE(Boson, X, 900001, 2, RETRIVE_PARAMETER(MX), nullptr);
REGISTER_PARTICLE(Boson, S, 900011, 1, nullptr, nullptr);

// * Particles entering the Boltzmann Equation
REGISTER_POI(900001, 1);

// * Regsiter Processes
REGISTER_PROCESS(XX_SS_Amp);

#endif  //_TOY_DM_H_s
