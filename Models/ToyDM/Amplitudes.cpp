#include "Amplitudes.h"

#include <cmath>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;

XX_HH_Amp::XX_HH_Amp() : Amplitude_Base("XX_SS") {
    Particle_Base *p_dm = RETRIEVE_PARTICLE(900001);
    Particle_Base *p_H = RETRIEVE_PARTICLE(900025);
    INITIAL.push_back(p_dm);
    INITIAL.push_back(p_dm);
    FINAL.push_back(p_H);
    FINAL.push_back(p_H);
    N_INITIAL = INITIAL.size();
    N_FINAL = FINAL.size();

    auto *ptr_lam = RETRIEVE_PARAMETER(Lam);
    Register_Dependencies(ptr_lam);
}

void XX_HH_Amp::Update_Value(REAL input) {
    REAL lam = RETRIEVE_PARAMETER(Lam)->Get_Value();
    Sub1 = pow(lam, 2);
}

void XX_HH_Amp::Update_Amp(REAL sqrt_shat) {
    REAL mH = RETRIEVE_PARTICLE(900025)->Get_Mass();
    REAL MH2 = mH * mH;
    REAL s = sqrt_shat * sqrt_shat;
    REAL kallen_sqrt = sqrt(Kallen_Lam(1.0, MH2 / s, MH2 / s));
    amp_res = Sub1 * kallen_sqrt / 8.0 / M_PI;
}

REAL XX_HH_Amp::Get_Coeff(REAL T, int PID) {
    if (PID == 900001) {
        Particle_Base *pp = Particle_Factory::Get_Particle_Factory().Get_Particle(900001);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        REAL res = 2.0 * (1.0 - pow(Y / YeqT, 2));
        if (fabs(res) < 1e-5) {
            res = 2.0 * (2 * pp->Delta_Yield_Ratio);
        }

        return res;
    } else {
        return 0;
    }
}
