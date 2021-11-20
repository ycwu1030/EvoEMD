#include "Amplitudes.h"
using namespace EvoEMD;
N_LPhi_Amp::N_LPhi_Amp() : Amplitude_Base("Decay_NLPhi") {
    Particle_Base *p_N1 = RETRIEVE_PARTICLE(900001);
    Particle_Base *p_l = RETRIEVE_PARTICLE(900011);
    Particle_Base *p_phi = RETRIEVE_PARTICLE(25);

    FINAL.push_back(p_l);
    FINAL.push_back(p_phi);
    INITIAL.push_back(p_N1);

    N_INITIAL = INITIAL.size();
    N_FINAL = FINAL.size();

    Parameter_Base *ptr_lam = RETRIEVE_PARAMETER(Lam);
    Parameter_Base *ptr_mn1 = RETRIEVE_PARAMETER(MN1);
    Register_Dependencies(ptr_lam, ptr_mn1);
}
void N_LPhi_Amp::Update_Value(REAL input) {
    REAL lam = GET_PARAMETER_VALUE(Lam);
    REAL mn1 = GET_PARAMETER_VALUE(MN1);

    Sub1 = 4 * lam * lam * mn1 * mn1;
}

void N_LPhi_Amp::Update_Amp(REAL sqrt_shat) { amp_res = 4 * M_PI * Sub1; }

REAL N_LPhi_Amp::Get_Coeff(REAL T, int PID) {
    if (PID == 900001) {
        Particle_Base *pp = RETRIEVE_PARTICLE(900001);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        REAL res = (1.0 - Y / YeqT);
        if (fabs(res) < 1e-5) {
            res = pp->Delta_Yield_Ratio;
        }
        return res;
    } else if (PID == 900011) {
        Particle_Base *pp = RETRIEVE_PARTICLE(900011);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        REAL res = -Y / YeqT / 2.0;
        return res;
    } else {
        return 0;
    }
}

delta_N_LPhi_Amp::delta_N_LPhi_Amp() : Amplitude_Base("Decay_decay_NLPhi") {
    Particle_Base *p_N1 = RETRIEVE_PARTICLE(900001);
    Particle_Base *p_l = RETRIEVE_PARTICLE(900011);
    Particle_Base *p_phi = RETRIEVE_PARTICLE(25);

    FINAL.push_back(p_l);
    FINAL.push_back(p_phi);
    INITIAL.push_back(p_N1);

    N_INITIAL = INITIAL.size();
    N_FINAL = FINAL.size();

    auto *ptr_lam = RETRIEVE_PARAMETER(Lam);
    auto *ptr_mn1 = RETRIEVE_PARAMETER(MN1);
    auto *ptr_eps = RETRIEVE_PARAMETER(Eps);
    Register_Dependencies(ptr_lam, ptr_mn1, ptr_eps);
}

void delta_N_LPhi_Amp::Update_Value(REAL input) {
    REAL lam = GET_PARAMETER_VALUE(Lam);
    REAL mn1 = GET_PARAMETER_VALUE(MN1);
    REAL eps = GET_PARAMETER_VALUE(Eps);

    Sub1 = 4 * eps * lam * lam * mn1 * mn1;
}

void delta_N_LPhi_Amp::Update_Amp(REAL sqrt_shat) { amp_res = 4 * M_PI * Sub1; }

REAL delta_N_LPhi_Amp::Get_Coeff(REAL T, int PID) {
    REAL res = 0;
    if (PID == 900001) {
        Particle_Base *pp = RETRIEVE_PARTICLE(900011);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        res = -Y / YeqT / 2.0;
    } else if (PID == 900011) {
        Particle_Base *pp = RETRIEVE_PARTICLE(900001);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        res = -(1.0 - Y / YeqT);
        if (fabs(res) < 1e-5) {
            res = -pp->Delta_Yield_Ratio;
        }
    }
    return res;
}
