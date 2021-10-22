#ifndef __PHASE_SPACE_H__
#define __PHASE_SPACE_H__

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

/**
 * @brief  The Kallen lambda function
 * @note
 * @param  x: REAL parameter
 * @param  y: REAL parameter
 * @param  z: REAL parameter
 * @retval x^2+y^2+z^2-2xy-2yz-2zx
 */
REAL Kallen_Lam(REAL x, REAL y, REAL z);

/**
 * @brief  Calculate the energy of a particle in a 1->2 system
 * @note   Assume a 1->2 system: a->i j, the invariant mass of a is sqrt(s), mi
 * and mj is the invariant mass of the other two particle, calculate the energy
 * of i in the rest frame of a
 * @param  sqrt_s: The invariant mass of the mother system
 * @param  mi: The mass of particle i
 * @param  mj: The mass of particle j
 * @retval The energy of particle i in COM frame of i+j
 */
REAL Ei(REAL sqrt_s, REAL mi, REAL mj);

/**
 * @brief  Calculate the magnitude of momentum of a particle in a 1->2 system
 * @note   Assume a 1->2 system: a->i j, the invariant mass of a is sqrt(s), mi
 * and mj is the invariant mass of the other two particle, calculate the energy
 * of either i or j in the rest frame of a
 * @param  sqrt_s: The invariant mass of the mother system
 * @param  mi: The mass of particle i
 * @param  mj: The mass of particle j
 * @retval The magnitude of momentum of particle i or j in COM frame of i+j
 */
REAL Pi(REAL sqrt_s, REAL mi, REAL mj);

}  // namespace EvoEMD

#endif  //__PHASE_SPACE_H__
