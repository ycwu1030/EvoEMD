#ifndef _TOY_DM_PARTICLES_H_
#define _TOY_DM_PARTICLES_H_

#include "EvoEMD/ParticleBase.h"
#include "Parameters.h"
using namespace EvoEMD;
REGISTER_PARTICLE(Boson, X, 900001, 2, RETRIVE_PARAMETER(MX), nullptr, true);
REGISTER_POI(900001, 1);

REGISTER_PARTICLE(Boson, S, 900011, 1, nullptr, nullptr, true);

#endif  //_TOY_DM_PARTICLES_H_
