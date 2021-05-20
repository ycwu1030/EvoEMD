#ifndef __LEPTOGENESIS_RATE_H__
#define __LEPTOGENESIS_RATE_H__

#include "EMD.h"
#include "Neutrino.h"

class LeptogenesisRate
{
public:
    LeptogenesisRate();
    ~LeptogenesisRate(){};

    void Set_Nu_Masses(REAL M1, REAL M2, REAL M3);
    void Set_Nu_MassOrdering(Nu_TypeI_SeeSaw::MassOrdering od = Nu_TypeI_SeeSaw::NORMAL_ORDER);
    void Set_Dark_Sector_Masses(REAL mchi, REAL ms);
    void Set_LambdaX(REAL lamx);

    void Update();

    int Get_DOF() { return DOF; }

    VD Ri(REAL Temp, VD y);

    REAL Get_Temperature_from_Z(REAL z)
    {
        return MNR1 / z;
    }
    REAL Get_Z_from_Temperature(REAL T)
    {
        return MNR1 / T;
    }

private:
    int DOF;
    VD Yeq;
    REAL MNR1, MNR2, MNR3;
    REAL MCHI, MS;
    REAL LamX;

    REAL GammaN1, GammaN2;
    REAL eps1, eps2;

    Nu_TypeI_SeeSaw Nu_Param;

    REAL gamma_Ni_Chi_S[2];
    REAL gamma_S_Chi_Ni[2];
    REAL gamma_Chi_S_Ni[2];
    REAL gamma_L_Phi_Chi_S;
    REAL gamma_Ni_Ni_Chi_Chi[2];
    REAL gamma_Ni_Ni_S_S[2];
    REAL gamma_Ni_L_Phi[2];
    void Calc_N_Width();
    void Calc_CP_Asym();
    void Calc_Gammas(REAL Temp);
    REAL Calc_NChiS_Gamma(REAL Temp, REAL MA, REAL MB, REAL MC); // A -> B C, for N Chi S system
    REAL Calc_NLPhi_Gamma(REAL Temp, int i);
    REAL Calc_LPhiChiS_Gamma(REAL Temp);
    REAL Calc_NNChiChi_Gamma(REAL Temp, int i);
    REAL Calc_NNSS_Gamma(REAL Temp, int i);

    void Calc_Yield_Equilibrium(REAL Temp);

    REAL R0(REAL Temp, VD y); // For Lepton Number;
    REAL R1(REAL Temp, VD y); // For N1
    REAL R2(REAL Temp, VD y); // For N2
    REAL R3(REAL Temp, VD y); // For Chi
    REAL R4(REAL Temp, VD y); // For S
};

REAL SQkaller(REAL x, REAL y, REAL z);

#endif //__LEPTOGENESIS_RATE_H__
