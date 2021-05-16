#ifndef EMD_H
#define EMD_H

#include "RealTypes.h"

class EMD
{
private:
    REAL TInflation; // Temperature at the end of Inflation.
    REAL Ti;         // Temperature at the beginning of Early Matter Dominant Era
    REAL Te;         // Temperature at which we start the Entropy Production period
    REAL Tr;         // Reheating Temperature at the end of the Early Matter Dominant Era (at the end of Entropy production)
    REAL Tf;         // Final Temperature where we stop the evolution.

    // D.O.F
    REAL gei, gsi;
    REAL ger, gsr;
    REAL gee, gse;

    // Parameter related to EMD:
    REAL Delta;
    REAL CoverD;
    REAL Hubble_RD_at_Tr;

    void Get_Te();
    void Set_Temperature(REAL _Ti, REAL _Tr, REAL _Tf = 1e-3, REAL _TInflation = 1e15);

public:
    EMD();
    EMD(REAL _Ti, REAL _Tr, REAL _Tf = 1e-3, REAL _TInflation = 1e15);
    ~EMD() = default;

    REAL Get_Hubble_at_T(REAL Temp);
    friend double Equation_For_LogTe(double logTe, void *param);
};

/**
 * @brief  Calculate the particle number density at equilibrium
 * @note
 * @param  T: The temperature
 * @param  M: The mass of the particle
 * @param  g: The internal dof of the particle
 * @retval The particle number density at equilibrium
 */
REAL Number_Density_Eq(REAL T, REAL M, REAL g);

/**
 * @brief  Calculate the Entropy density at Temperature T
 * @note
 * @param  T: The Temperature
 * @retval The entropy density
 */
REAL Entropy_Density(REAL T);

/**
 * @brief Calculate the Yield (n/s) at equilibrium
 * @note
 * @param  T: The Temperature
 * @param  M: The mass of the particle
 * @param  g: The internal dof of the particle
 * @retval The Yield at equilibrium
 */
REAL Yield_Eq(REAL T, REAL M, REAL g);

#endif
