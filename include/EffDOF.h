#ifndef EFFDOF_H
#define EFFDOF_H

#include "RealTypes.h"

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

#endif //EFFDOF_H
