#include "Amplitudes.h"
using namespace EvoEMD;
N_LPhi_Amp::N_LPhi_Amp() : Amplitude() {
    Pseudo_Particle *p_N1 = RETRIVE_PARTICLE(900001);
    Pseudo_Particle *p_l = RETRIVE_PARTICLE(900011);
    Pseudo_Particle *p_phi = RETRIVE_PARTICLE(25);

    FINAL.push_back(p_l);
    FINAL.push_back(p_phi);
    INITIAL.push_back(p_N1);

    N_INITIAL = 1;
    N_FINAL = 2;

    amp_res.n_diag = 1;
    Process_Amp::CTH_RES_FULL den(1, 1);
    Process_Amp::PROPAGATOR_STRUCTURE ds = std::make_pair(0, den);
    amp_res.amps_denominator.push_back(std::make_pair(ds, ds));
    Process_Amp::NUMERATOR_STRUCTURE num(3, 0);
    amp_res.amps_numerator.push_back(num);
}

void N_LPhi_Amp::Update_Amp(REAL sqrt_shat) {
    REAL lam = GET_PARAMETER_VALUE(Lam);
    REAL mn1 = GET_PARAMETER_VALUE(MN1);
    amp_res.amps_numerator[0][0] = 4 * lam * lam * mn1 * mn1;
}

REAL N_LPhi_Amp::Get_Coeff(REAL T, int PID) {
    if (PID == 900001) {
        Pseudo_Particle *pp = RETRIVE_PARTICLE(900001);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        REAL res = (1.0 - Y / YeqT);
        if (fabs(res) < 1e-5) {
            res = pp->Delta_Yield_Ratio;
        }
        return res;
    } else if (PID == 900011) {
        Pseudo_Particle *pp = RETRIVE_PARTICLE(900011);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        REAL res = -Y / YeqT / 2.0;
        return res;
    } else {
        return 0;
    }
}

delta_N_LPhi_Amp::delta_N_LPhi_Amp() : Amplitude() {
    Pseudo_Particle *p_N1 = RETRIVE_PARTICLE(900001);
    Pseudo_Particle *p_l = RETRIVE_PARTICLE(900011);
    Pseudo_Particle *p_phi = RETRIVE_PARTICLE(25);

    FINAL.push_back(p_l);
    FINAL.push_back(p_phi);
    INITIAL.push_back(p_N1);

    N_INITIAL = 1;
    N_FINAL = 2;

    amp_res.n_diag = 1;
    Process_Amp::CTH_RES_FULL den(1, 1);
    Process_Amp::PROPAGATOR_STRUCTURE ds = std::make_pair(0, den);
    amp_res.amps_denominator.push_back(std::make_pair(ds, ds));
    Process_Amp::NUMERATOR_STRUCTURE num(3, 0);
    amp_res.amps_numerator.push_back(num);
}

void delta_N_LPhi_Amp::Update_Amp(REAL sqrt_shat) {
    REAL lam = GET_PARAMETER_VALUE(Lam);
    REAL mn1 = GET_PARAMETER_VALUE(MN1);
    REAL eps = GET_PARAMETER_VALUE(Eps);
    amp_res.amps_numerator[0][0] = 4 * eps * lam * lam * mn1 * mn1;
}

REAL delta_N_LPhi_Amp::Get_Coeff(REAL T, int PID) {
    REAL res = 0;
    if (PID == 900001) {
        Pseudo_Particle *pp = RETRIVE_PARTICLE(900011);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        res = -Y / YeqT / 2.0;
    } else if (PID == 900011) {
        Pseudo_Particle *pp = RETRIVE_PARTICLE(900001);
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
