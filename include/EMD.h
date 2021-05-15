#ifndef EMD_H
#define EMD_H

#include "RealTypes.h"

class EMD
{
private:
    REAL TInflation; // Temperature at the end of Inflation.
    REAL Ti; // Temperature at the beginning of Early Matter Dominant Era 
    REAL Te; // Temperature at which we start the Entropy Production period
    REAL Tr; // Reheating Temperature at the end of the Early Matter Dominant Era (at the end of Entropy production)
    REAL Tf; // Final Temperature where we stop the evolution.

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


/*
 * @brief: Calculate the number density at equilibrium
 * @param T, the temperature
 * @param M, the mass of the species
 * @param g, the internal dof of the species 
*/
REAL Number_Density_Eq(REAL T, REAL M, REAL g);

REAL Entropy_Density(REAL T);

REAL Yield_Eq(REAL T, REAL M, REAL g);

#endif