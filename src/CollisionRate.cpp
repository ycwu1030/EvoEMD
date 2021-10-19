#include "CollisionRate.h"

#include <cmath>

#include "PhaseSpace.h"
#include "gsl/gsl_sf_bessel.h"

REAL Decay12_Rate::Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) {
    // * For  1 -> 2 3;
    REAL sqrtshat = sqrt_shat;  // avoid unused warning
    REAL m1 = amp->INITIAL[0]->Get_Mass();
    REAL m2 = amp->FINAL[0]->Get_Mass();
    REAL m3 = amp->FINAL[1]->Get_Mass();
    REAL m12 = m1 * m1;
    REAL m22 = m2 * m2;
    REAL m32 = m3 * m3;

    // * Final state momentum when 1 is at rest;
    REAL pf = sqrt(Kallen_Lam(m12, m22, m32)) / 2.0 / m1;

    // * Two body phase space volumn
    REAL PS2 = 1.0 / (2.0 * m1) * pf / (16.0 * M_PI * M_PI * m1);
    REAL RES = 2.0                                     // * From integral over theta
               * 2.0 * M_PI                            // * From integral over phi
               * (amp->Get_Amp().Amp_Numerator[0][0])  // * The amplitude
               * PS2;                                  // * The two body phase space
    return RES;
}

REAL Decay12_Rate::Get_Collision_Rate(REAL T) {
    REAL AMP_WITH_PS = Get_Amp_Integrate_over_Phase_Space(0);

    REAL m1 = amp->INITIAL[0]->Get_Mass();
    REAL z = m1 / T;
    REAL k1z = gsl_sf_bessel_Kn(1, z);
    REAL RES = m1 * m1 * T / 2.0 / M_PI / M_PI * k1z;
    return RES;
}
