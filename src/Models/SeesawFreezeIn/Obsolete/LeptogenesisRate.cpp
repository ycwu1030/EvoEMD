#include "EvoEMD/LeptogenesisRate.h"

#include "EvoEMD/PhaseSpace.h"
#include "EvoEMD/spdlog_wrapper.h"

namespace EvoEMD {
LeptogenesisRate::LeptogenesisRate()
    : DOF(5), MNR1(3000), MNR2(3000), MNR3(1e10), MCHI(100), MS(300), LamX(1.0), Yeq(5) {
    Nu_Param.Set_Heavy_Neutrino_Mass(MNR1, MNR2, MNR3);
}

void LeptogenesisRate::Set_Nu_Masses(REAL M1, REAL M2, REAL M3) {
    MNR1 = M1;
    MNR2 = M2;
    MNR3 = M3;
    Nu_Param.Set_Heavy_Neutrino_Mass(MNR1, MNR2, MNR3);
}

void LeptogenesisRate::Set_Nu_MassOrdering(Nu_TypeI_SeeSaw::MassOrdering od) { Nu_Param.Set_Mass_Ordering(od); }

void LeptogenesisRate::Set_Dark_Sector_Masses(REAL mchi, REAL ms) {
    MCHI = mchi;
    MS = ms;
}

void LeptogenesisRate::Set_LambdaX(REAL lamx) { LamX = lamx; }

void LeptogenesisRate::Calc_N_Width() {
    GammaN1 = MNR1 / 8 / M_PI * real(Nu_Param.Get_YdagYij(1, 1));
    GammaN2 = MNR2 / 8 / M_PI * real(Nu_Param.Get_YdagYij(2, 2));
    if (MNR1 > MCHI + MS) {
        GammaN1 += MNR1 / 32.0 / M_PI * 2 * LamX * LamX * (1.0 + (MCHI * MCHI - MS * MS) / MNR1 / MNR1) *
                   sqrt(Kallen_Lam(1.0, MCHI * MCHI / MNR1 / MNR1, MS * MS / MNR1 / MNR1));
    }
    if (MNR2 > MCHI + MS) {
        GammaN2 += MNR2 / 32.0 / M_PI * 2 * LamX * LamX * (1.0 + (MCHI * MCHI - MS * MS) / MNR2 / MNR2) *
                   sqrt(Kallen_Lam(1.0, MCHI * MCHI / MNR2 / MNR2, MS * MS / MNR2 / MNR2));
    }
    SPDLOG_INFO_FILE("Width for N1: GammaN1 = {:+9.8e}.", GammaN1);
    SPDLOG_INFO_FILE("Width for N2: GammaN2 = {:+9.8e}.", GammaN2);
}

void LeptogenesisRate::Calc_CP_Asym() {
    auto yy12 = Nu_Param.Get_YdagYij(1, 2);
    auto yy21 = Nu_Param.Get_YdagYij(2, 1);
    REAL Iyy122 = imag(yy12 * yy12);
    REAL Iyy212 = imag(yy21 * yy21);
    REAL yy11 = real(Nu_Param.Get_YdagYij(1, 1));
    REAL yy22 = real(Nu_Param.Get_YdagYij(2, 2));
    REAL DM2 = MNR1 * MNR1 - MNR2 * MNR2;
    eps1 = Iyy122 / yy11 / yy22 * DM2 * MNR1 * GammaN2 / (pow(DM2, 2) + pow(MNR1 * GammaN2, 2));
    eps2 = Iyy212 / yy22 / yy11 * (-DM2) * MNR2 * GammaN1 / (pow(DM2, 2) + pow(MNR2 * GammaN1, 2));
    SPDLOG_INFO_FILE("CP asymmetry for N1: epsilon1 = {:+9.8e}.", eps1);
    SPDLOG_INFO_FILE("CP asymmetry for N2: epsilon2 = {:+9.8e}.", eps2);
}

void LeptogenesisRate::Update() {
    SPDLOG_INFO_FILE("Updating Width and CP asymmetry.");
    Calc_N_Width();
    Calc_CP_Asym();
}

VD LeptogenesisRate::Ri(REAL Temp, VD y) {
    Calc_Yield_Equilibrium(Temp);
    Calc_Gammas(Temp);
    SPDLOG_INFO_FILE("Calculate Collision Rate at T = {:+9.8e}", Temp);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_INFO_FILE("                         at Y[{}] = {:+9.8e}", i, y[i]);
    }
    return {R0(Temp, y), R1(Temp, y), R2(Temp, y), R3(Temp, y), R4(Temp, y)};
}

void LeptogenesisRate::Calc_Gammas(REAL Temp) {
    gamma_Ni_Chi_S[0] = Calc_NChiS_Gamma(Temp, 0);
    gamma_Ni_Chi_S[1] = Calc_NChiS_Gamma(Temp, 1);
    gamma_S_Chi_Ni[0] = Calc_SChiN_Gamma(Temp, 0);
    gamma_S_Chi_Ni[1] = Calc_SChiN_Gamma(Temp, 1);
    gamma_Chi_S_Ni[0] = Calc_ChiSN_Gamma(Temp, 0);
    gamma_Chi_S_Ni[1] = Calc_ChiSN_Gamma(Temp, 1);
    gamma_L_Phi_Chi_S = Calc_LPhiChiS_Gamma(Temp);
    gamma_Ni_Ni_Chi_Chi[0][0] = Calc_NNChiChi_Gamma(Temp, 0, 0);
    gamma_Ni_Ni_Chi_Chi[0][1] = Calc_NNChiChi_Gamma(Temp, 0, 1);
    gamma_Ni_Ni_Chi_Chi[1][0] = Calc_NNChiChi_Gamma(Temp, 1, 0);
    gamma_Ni_Ni_Chi_Chi[1][1] = Calc_NNChiChi_Gamma(Temp, 1, 1);
    gamma_Ni_Ni_S_S[0][0] = Calc_NNSS_Gamma(Temp, 0, 0);
    gamma_Ni_Ni_S_S[0][1] = Calc_NNSS_Gamma(Temp, 0, 1);
    gamma_Ni_Ni_S_S[1][0] = Calc_NNSS_Gamma(Temp, 1, 0);
    gamma_Ni_Ni_S_S[1][1] = Calc_NNSS_Gamma(Temp, 1, 1);
    gamma_Ni_L_Phi[0] = Calc_NLPhi_Gamma(Temp, 0);
    gamma_Ni_L_Phi[1] = Calc_NLPhi_Gamma(Temp, 1);
}

void LeptogenesisRate::Calc_Yield_Equilibrium(REAL Temp) {
    Yeq[0] = Yield_Eq_FD(Temp, 2 * 2);     // For lepton
    Yeq[1] = Yield_Eq(Temp, MNR1, 2);      // For N1;
    Yeq[2] = Yield_Eq(Temp, MNR2, 2);      // For N2;
    Yeq[3] = Yield_Eq(Temp, MCHI, 2 * 2);  // For Chi
    Yeq[4] = Yield_Eq(Temp, MS, 2);        // For S
    SPDLOG_INFO_FILE("Calculate Yield at Equilibrium at T = {:+9.8}.", Temp);
    SPDLOG_INFO_FILE("For Lep, Yeq = {:+9.8e}", Yeq[0]);
    SPDLOG_INFO_FILE("For N1,  Yeq = {:+9.8e}", Yeq[1]);
    SPDLOG_INFO_FILE("For N2,  Yeq = {:+9.8e}", Yeq[2]);
    SPDLOG_INFO_FILE("For Chi, Yeq = {:+9.8e}", Yeq[3]);
    SPDLOG_INFO_FILE("For S,   Yeq = {:+9.8e}", Yeq[4]);
}

REAL LeptogenesisRate::R0(REAL Temp, VD y) {
    REAL R0 = eps1 * gamma_Ni_L_Phi[0] * (y[1] / Yeq[1] - 1.0) - gamma_Ni_L_Phi[0] * y[0] / 2.0 / Yeq[0];
    R0 += eps2 * gamma_Ni_L_Phi[1] * (y[2] / Yeq[2] - 1.0) - gamma_Ni_L_Phi[1] * y[0] / 2.0 / Yeq[0];
    // TODO: Adding terms related to 2 -> 2 processes
    SPDLOG_INFO_FILE("Collision Rate for Lepton Number: R0 = {:+9.8e}", R0);
    return R0;
}

REAL LeptogenesisRate::R1(REAL Temp, VD y) {
    REAL RN1 = -gamma_Ni_L_Phi[0] * (y[1] / Yeq[1] - 1.0);
    RN1 += -eps1 * gamma_Ni_L_Phi[0] * y[0] / 2 / Yeq[0];
    RN1 += -gamma_Ni_Ni_Chi_Chi[0][0] * (y[1] * y[1] / Yeq[1] / Yeq[1] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RN1 += -gamma_Ni_Ni_Chi_Chi[0][1] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RN1 += -gamma_Ni_Ni_Chi_Chi[1][0] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RN1 += -gamma_Ni_Ni_S_S[0][0] * (y[1] * y[1] / Yeq[1] / Yeq[1] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RN1 += -gamma_Ni_Ni_S_S[0][1] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RN1 += -gamma_Ni_Ni_S_S[1][0] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    // if (MNR1 > MCHI + MS)
    // { // * The mass threshold will be checked when calculate gamma
    RN1 += -gamma_Ni_Chi_S[0] * (y[1] / Yeq[1] - y[3] * y[4] / Yeq[3] / Yeq[4]);
    // }
    SPDLOG_INFO_FILE("Collision Rate for N1: R1 = {:+9.8e}", RN1);
    return RN1;
}

REAL LeptogenesisRate::R2(REAL Temp, VD y) {
    REAL RN2 = -gamma_Ni_L_Phi[1] * (y[2] / Yeq[2] - 1.0);
    RN2 += -eps2 * gamma_Ni_L_Phi[1] * y[0] / 2 / Yeq[0];
    RN2 += -gamma_Ni_Ni_Chi_Chi[1][1] * (y[2] * y[2] / Yeq[2] / Yeq[2] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RN2 += -gamma_Ni_Ni_Chi_Chi[1][0] * (y[2] * y[1] / Yeq[2] / Yeq[1] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RN2 += -gamma_Ni_Ni_Chi_Chi[0][1] * (y[2] * y[1] / Yeq[2] / Yeq[1] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RN2 += -gamma_Ni_Ni_S_S[1][1] * (y[2] * y[2] / Yeq[2] / Yeq[2] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RN2 += -gamma_Ni_Ni_S_S[1][0] * (y[2] * y[1] / Yeq[2] / Yeq[1] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RN2 += -gamma_Ni_Ni_S_S[0][1] * (y[2] * y[1] / Yeq[2] / Yeq[1] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    // if (MNR2 > MCHI + MS)
    // {
    RN2 += -gamma_Ni_Chi_S[1] * (y[2] / Yeq[2] - y[3] * y[4] / Yeq[3] / Yeq[4]);
    // }
    SPDLOG_INFO_FILE("Collision Rate for N2: R2 = {:+9.8e}", RN2);
    return RN2;
}

REAL LeptogenesisRate::R3(REAL Temp, VD y) {
    REAL RCHI = gamma_L_Phi_Chi_S;
    RCHI += gamma_Ni_Ni_Chi_Chi[0][0] * (y[1] * y[1] / Yeq[1] / Yeq[1] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RCHI += gamma_Ni_Ni_Chi_Chi[0][1] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RCHI += gamma_Ni_Ni_Chi_Chi[1][0] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    RCHI += gamma_Ni_Ni_Chi_Chi[1][1] * (y[2] * y[2] / Yeq[2] / Yeq[2] - y[3] * y[3] / Yeq[3] / Yeq[3]);
    // if (MNR1 > MCHI + MS)
    // {
    RCHI += gamma_Ni_Chi_S[0] * (y[1] / Yeq[1] - y[3] * y[4] / Yeq[3] / Yeq[4]);
    // }
    // if (MNR2 > MCHI + MS)
    // {
    RCHI += gamma_Ni_Chi_S[1] * (y[2] / Yeq[2] - y[3] * y[4] / Yeq[3] / Yeq[4]);
    // }
    RCHI += gamma_S_Chi_Ni[0] * (y[4] / Yeq[4] - y[3] * y[1] / Yeq[3] / Yeq[1]);
    RCHI += gamma_S_Chi_Ni[1] * (y[4] / Yeq[4] - y[3] * y[2] / Yeq[3] / Yeq[2]);
    RCHI += -gamma_Chi_S_Ni[0] * (y[3] / Yeq[3] - y[4] * y[1] / Yeq[4] / Yeq[1]);
    RCHI += -gamma_Chi_S_Ni[1] * (y[3] / Yeq[3] - y[4] * y[2] / Yeq[4] / Yeq[2]);
    SPDLOG_INFO_FILE("Collision Rate for CHI: R3 = {:+9.8e}", RCHI);
    return RCHI;
}

REAL LeptogenesisRate::R4(REAL Temp, VD y) {
    REAL RS = gamma_L_Phi_Chi_S;
    RS += gamma_Ni_Ni_S_S[0][0] * (y[1] * y[1] / Yeq[1] / Yeq[1] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RS += gamma_Ni_Ni_S_S[0][1] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RS += gamma_Ni_Ni_S_S[1][0] * (y[1] * y[2] / Yeq[1] / Yeq[2] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RS += gamma_Ni_Ni_S_S[1][1] * (y[2] * y[2] / Yeq[2] / Yeq[2] - y[4] * y[4] / Yeq[4] / Yeq[4]);
    RS += gamma_Ni_Chi_S[0] * (y[1] / Yeq[1] - y[3] * y[4] / Yeq[3] / Yeq[4]);
    RS += gamma_Ni_Chi_S[1] * (y[2] / Yeq[2] - y[3] * y[4] / Yeq[3] / Yeq[4]);
    RS += -gamma_S_Chi_Ni[0] * (y[4] / Yeq[4] - y[3] * y[1] / Yeq[3] / Yeq[1]);
    RS += -gamma_S_Chi_Ni[1] * (y[4] / Yeq[4] - y[3] * y[2] / Yeq[3] / Yeq[2]);
    RS += gamma_Chi_S_Ni[0] * (y[3] / Yeq[3] - y[4] * y[1] / Yeq[4] / Yeq[1]);
    RS += gamma_Chi_S_Ni[1] * (y[3] / Yeq[3] - y[4] * y[2] / Yeq[4] / Yeq[2]);
    SPDLOG_INFO_FILE("Collision Rate for S: R4 = {:+9.8e}", RS);
    return RS;
}
}  // namespace EvoEMD
