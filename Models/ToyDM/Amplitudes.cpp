#include "Amplitudes.h"

#include <cmath>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;

XX_SS_Amp::XX_SS_Amp() : Amplitude_Base("XX_SS") {
    Particle_Base *p_dm = RETRIEVE_PARTICLE(900001);
    Particle_Base *p_s = RETRIEVE_PARTICLE(900011);
    INITIAL.push_back(p_dm);
    INITIAL.push_back(p_dm);
    FINAL.push_back(p_s);
    FINAL.push_back(p_s);
    N_INITIAL = INITIAL.size();
    N_FINAL = FINAL.size();

    auto *ptr_lam = RETRIEVE_PARAMETER(Lam);
    Register_Dependencies(ptr_lam);
}

void XX_SS_Amp::Update_Value(REAL input) {
    REAL lam = RETRIEVE_PARAMETER(Lam)->Get_Value();
    Sub1 = pow(lam, 4);
}

void XX_SS_Amp::Update_Amp(REAL sqrt_shat) {
    // * final states s is massless, so sqrt(lam(1,0,0))=1
    amp_res = Sub1 / 8.0 / M_PI;
}

REAL XX_SS_Amp::Get_Coeff(REAL T, int PID) {
    if (PID == 900001) {
        Particle_Base *pp = Particle_Factory::Get_Particle_Factory().Get_Particle(900001);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        REAL res = (1.0 - pow(Y / YeqT, 2));
        if (fabs(res) < 1e-5) {
            res = 2 * pp->Delta_Yield_Ratio;
        }

        return res;
    } else {
        return 0;
    }
}
