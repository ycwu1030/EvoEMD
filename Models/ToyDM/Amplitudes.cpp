#include "Amplitudes.h"

#include <cmath>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;

XX_SS_Amp::XX_SS_Amp() : Amplitude(1) {
    Particle_Factory &pf = Particle_Factory::Get_Particle_Factory();
    Pseudo_Particle *p_dm = pf.Get_Particle(900001);
    Pseudo_Particle *p_l = pf.Get_Particle(900011);
    FINAL.push_back(p_dm);
    FINAL.push_back(p_dm);
    INITIAL.push_back(p_l);
    INITIAL.push_back(p_l);
    N_INITIAL = 2;
    N_FINAL = 2;
    // * For the diagram-0
    auto &amp_0 = amp_res[0];
    auto ds = Build_Propagator(0, 1.0);
    amp_0.Denominator = Build_Denominator(ds, ds);
    amp_0.Numerator = Build_Numerator(1.0, 0.0, 0.0);
}

void XX_SS_Amp::Update_Amp(REAL sqrt_shat) {
    REAL lam = RETRIEVE_PARAMETER(Lam)->Get_Value();

    // * Update amp for diagram-0
    amp_res[0].Numerator = Build_Numerator(lam * lam * lam * lam, 0.0, 0.0);
}

REAL XX_SS_Amp::Get_Coeff(REAL T, int PID) {
    if (PID == 900001) {
        Pseudo_Particle *pp = Particle_Factory::Get_Particle_Factory().Get_Particle(900001);
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
