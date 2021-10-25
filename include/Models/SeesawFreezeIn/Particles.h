#ifndef _SEESAW_FI_PARTICLES_H_
#define _SEESAW_FI_PARTICLES_H_

#include "EvoEMD/ParticleBase.h"
#include "Parameters.h"
using namespace EvoEMD;

REGISTER_PARTICLE(Fermion, N1, 900001, 2, RETRIVE_PARAMETER(MN1), nullptr, true);
REGISTER_POI(900001, 1);

REGISTER_PARTICLE(Fermion, N2, 900002, 2, RETRIVE_PARAMETER(MN2), nullptr, true);

REGISTER_PARTICLE(Fermion, N3, 900003, 2, RETRIVE_PARAMETER(MN3), nullptr, true);

#endif  //_SEESAW_FI_PARTICLES_H_
