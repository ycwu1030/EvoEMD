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
REAL f_ge(REAL T);

/**
 * @brief  Calculate the effective degree of freedom for entropy
 * @note   This is for SM
 * @param  T: The temperature
 * @retval The effective dof at T
 */
REAL f_gs(REAL T);

/**
 * @brief  Calculate ge_star = 1 + dlog(ge)/dlog(T) /4;
 * @note   This is for SM
 * @param  T: The temperature
 * @retval ge_star
 */
REAL f_ge_star(REAL T);

/**
 * @brief  Calculate gs_star = 1 + dlog(gs)/dlog(T) /3;
 * @note   This is for SM
 * @param  T: The temperature
 * @retval gs_star
 */
REAL f_gs_star(REAL T);

}  // namespace EvoEMD

#endif  // EFFDOF_H
