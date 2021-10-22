#ifndef _EFFDOF_H_
#define _EFFDOF_H_

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

/**
 * @brief  Calculate the effective degree of freedom for energy density
 * @note   This is for SM
 * @param  T: The temperature
 * @retval The effective dof at T
 */
REAL ge(REAL T);

/**
 * @brief  Calculate the effective degree of freedom for entropy
 * @note   This is for SM
 * @param  T: The temperature
 * @retval The effective dof at T
 */
REAL gs(REAL T);

}  // namespace EvoEMD

#endif  // EFFDOF_H
