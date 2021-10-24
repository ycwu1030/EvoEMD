#ifndef EMD_H
#define EMD_H

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

class EMD {
private:
    REAL TInflation;  // Temperature at the end of Inflation.
    REAL Ti;          // Temperature at the beginning of Early Matter Dominant Era
    REAL Te;          // Temperature at which we start the Entropy Production period
    REAL Tr;          // Reheating Temperature at the end of the Early Matter Dominant
                      // Era (at the end of Entropy production)
    REAL Tf;          // Final Temperature where we stop the evolution.

    // D.O.F
    REAL gei, gsi;
    REAL ger, gsr;
    REAL gee, gse;

    // Parameter related to EMD:
    REAL Delta;
    REAL CoverD;  // bla
    REAL Hubble_RD_at_Tr;

    /**
     * @brief  Solve the Te from Ti and Tr
     * @note   Only use internally
     * @retval None, setting Te internally
     */
    void Solve_Te();

public:
    typedef enum { ERDE = 0, EMDE = 1, EPE = 2, RDE = 3 } Period;
    EMD();

    /**
     * @brief Constructor with specified Ti, Tr, Tf and TInflation
     * @note
     * @param  _Ti: The initial temperature of EMD
     * @param  _Tr: The end temperature of the EMD/EPE
     * @param  _Tf: Default(1e-3), the smallest temperature that we will end any
     * calculation
     * @param  _TInflation: Default(1e15), the highest temperature that we will
     * consider
     * @retval
     */
    EMD(REAL _Ti, REAL _Tr, REAL _Tf = 1e-3, REAL _TInflation = 1e15);
    ~EMD() = default;

    /**
     * @brief  Set the corresponding temperature
     * @note
     * @param  _Ti: The initial temperature of EMD
     * @param  _Tr: The end temperature of the EMD/EPE
     * @param  _Tf: Default(1e-3), the smallest temperature that we will end any
     * calculation
     * @param  _TInflation: Default(1e15), the highest temperature that we will
     * consider
     * @retval None
     */
    void Set_Temperature(REAL _Ti, REAL _Tr, REAL _Tf = 1e-3, REAL _TInflation = 1e15);

    /**
     * @brief  Get the Hubble parameter at Temperature Temp
     * @note   All the calculation should be performed in Natural unit,
     *         and Temperature is in GeV
     *         We will choose the period according to Temp
     * @param  Temp: The temperature in GeV
     * @retval The Hubble parameter at Temp
     */
    REAL Get_Hubble_at_T(REAL Temp);

    /**
     * @brief  Get the Hubble parameter for Period pd at Temperature Temp
     * @note   All the calculation should be performed in Natural unit,
     *         and Temperature is in GeV.
     *         The Hubble parameter will always be calculated using the
     *         formula for pd, even if Temp is outside the corresponding pd
     * @param  pd: The period (ERDE, EMDE, EPE, RDE)
     * @param  Temp: The temperature in GeV
     * @retval
     */
    REAL Get_Hubble_at_T_Period(Period pd, REAL Temp);

    /**
     * @brief  Get the Delta, which indicates how long the EMD last.
     *         Delta also indicates how much entropy is produced.
     * @note
     * @retval The Delta value;
     */
    REAL Get_Delta() const { return Delta; }
    REAL Get_Ti() const { return Ti; }
    REAL Get_Tr() const { return Tr; }
    REAL Get_Te() const { return Te; }
    REAL Get_TInflation() const { return TInflation; }
    REAL Get_Tf() const { return Tf; }

    // * This is the function used to obtain the solution Te,
    // * shouldn't be usedexternally
    friend double Equation_For_LogTe(double logTe, void *param);
};

/**
 * @brief  Calculate the particle number density at equilibrium
 * @note   with _Massless is for massless case.
 *         Note that we use the Maxwell distribution here, not
 *         Boltzmann-Einstein or Fermi-Dirac distribution.
 * @param  T: The temperature
 * @param  M: The mass of the particle
 * @param  g: The internal dof of the particle
 * @retval The particle number density at equilibrium
 */
REAL Number_Density_Eq(REAL T, REAL M, REAL g);  // Equilibrium number density - Maxwell-Boltzmann statistics
REAL Number_Density_Eq_BE(REAL T, REAL g);       // Equilibrium number density - Bose-Einstein statistics
REAL Number_Density_Eq_FD(REAL T, REAL g);       // Equilibrium number density - Fermi-Dirac statistics

// TODO: Using the Boltzmann-Einstein/Fermi-Dirac distribution for the number
// density

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
REAL Yield_Eq_BE(REAL T, REAL g);
REAL Yield_Eq_FD(REAL T, REAL g);

}  // namespace EvoEMD

#endif
