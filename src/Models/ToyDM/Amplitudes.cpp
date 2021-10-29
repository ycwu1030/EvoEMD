#include "Models/ToyDM/Amplitudes.h"

#include <cmath>

#include "EvoEMD/ParticleBase.h"

using namespace EvoEMD;

XX_SS_Amp::XX_SS_Amp() : Amplitude() {
    Particle_Factory &pf = Particle_Factory::Get_Particle_Factory();
    Pseudo_Particle *p_dm = pf.Get_Particle(900001);
    Pseudo_Particle *p_l = pf.Get_Particle(900011);
    FINAL.push_back(p_dm);
    FINAL.push_back(p_dm);
    INITIAL.push_back(p_l);
    INITIAL.push_back(p_l);
    N_INITIAL = 2;
    N_FINAL = 2;
    amp_res.n_diag = 1;
    Process_Amp::CTH_RES_FULL den(1, 1);
    Process_Amp::PROPAGATOR_STRUCTURE ds = std::make_pair(0, den);
    amp_res.amps_denominator.push_back(std::make_pair(ds, ds));
    Process_Amp::NUMERATOR_STRUCTURE num(3, 0);
    amp_res.amps_numerator.push_back(num);
}

void XX_SS_Amp::Update_Amp(REAL sqrt_shat) {
    REAL lam = RETRIVE_PARAMETER(Lam)->Get_Value();
    amp_res.amps_numerator[0][0] = lam * lam * lam * lam;
}

REAL XX_SS_Amp::Get_Coeff(REAL T, int PID) {
    if (PID == 900001) {
        Pseudo_Particle *pp = Particle_Factory::Get_Particle_Factory().Get_Particle(900001);
        REAL Y = pp->Yield;
        REAL YeqT = pp->Get_Equilibrium_Yield_at_T(T);
        if (YeqT == 0) return 0;
        REAL res = (1.0 - pow(Y / YeqT, 2));
        if (res < 1e-5) {
            res = 2 * pp->Delta_Yield_Ratio;
        }

        return res;
    } else {
        return 0;
    }
}
