#include "CollisionRate.h"

#include <cmath>

#include "PhaseSpace.h"
#include "gsl/gsl_sf_bessel.h"

REAL Decay12_Rate::Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) {
    // * For  1 -> 2 3;
    REAL m1 = amp->INITIAL[0]->Get_Mass();
    REAL m2 = amp->FINAL[0]->Get_Mass();
    REAL m3 = amp->FINAL[1]->Get_Mass();
    REAL m12 = m1 * m1;
    REAL m22 = m2 * m2;
    REAL m32 = m3 * m3;

    // * Final state momentum when 1 is at rest;
    REAL pf = sqrt(Kallen_Lam(m12, m22, m32)) / 2.0 / m1;

    // * Two body phase space volumn
    const Process_Amp &res = amp->Get_Amp(sqrt_shat);
    REAL PS2 = 1.0 / (2.0 * m1) * pf / (16.0 * M_PI * M_PI * m1);
    REAL RES = 2.0                                // * From integral over cos(theta)
               * 2.0 * M_PI                       // * From integral over phi
               * (res.Get_Numerator()[0].second)  // * The amplitude
               * PS2;                             // * The two body phase space
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

REAL Scatter22_Rate::Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) {
    using NS = Process_Amp::NUMERATOR_STRUCTURE;
    using DS = Process_Amp::DENOMINATOR_STRUCTURE;
    const Process_Amp &amp_res = amp->Get_Amp(sqrt_shat);
    REAL amp_total = 0;
    for (auto &&chan : Process_Amp::All_Channel) {
        unsigned n_diag = amp_res.n_diag.at(chan);
        for (unsigned i_diag = 0; i_diag < n_diag; i_diag++) {
            NS &num = amp_res.Get_Numerator(chan, i_diag);
            DS &den = amp_res.Get_Denominator(chan, i_diag);
            amp_total += Get_Amp_Integrate_over_Phase_Space_Single_Channel(num, den);
        }
    }
    return amp_total;
}

inline REAL INT_RES_00(REAL a, REAL b, REAL c, REAL d, REAL e) {
    // * Int[(a+b x+ c x^2)/d/e,{x,-1,1}]
    return 2 * (3 * a + c) / 3.0 / d / e;
}
inline REAL INT_RES_01(REAL a, REAL b, REAL c, REAL d, REAL e, REAL f) {
    // * Int[(a+b x+ c x^2)/(d+ e x)/f, {x,-1,1}]
    REAL den = pow(e, 3) * f;
    REAL num = 2 * e * (b * e - c * d) - log((d - e) / (d + e)) * (e * (a * e - b * d) + c * d * d);
    return num / den;
}

inline REAL INT_RES_11_0(REAL a, REAL b, REAL c, REAL d, REAL e, REAL f, REAL g) {
    // * Int[(a+b x+c x^2)/(d+e x)/(f+g x),{x,-1,1}]
    REAL den = e * e * g * g * (e * f - d * g);
    REAL num1 = 2 * c * e * g * (e * f - d * g);
    REAL num2 = g * g * (a * e * e - b * d * e + c * d * d) * log((d + e) / (d - e));
    REAL num3 = e * e * (a * g * g - b * f * g + c * f * f) * log((f - g) / (f + g));
    return (num1 + num2 + num3) / den;
}

inline REAL INT_RES_11_1(REAL a, REAL b, REAL c, REAL d, REAL e) {
    // * Int[(a+b x+c x^2)/(d+ e x)^2,{x,-1,1}]
    REAL den = pow(e, 3) * (e - d) * (d + e);
    REAL num = 2 * e * (e * e * (c - a) + b * d * e - 2 * c * d * d) -
               (d - e) * (d + e) * log((d - e) / (d + e)) * (2 * c * d - b * e);
    return num / den;
}

REAL Scatter22_Rate::Get_Amp_Integrate_over_Phase_Space_Single_Channel(
    const Process_Amp::NUMERATOR_STRUCTURE &numerator, const Process_Amp::DENOMINATOR_STRUCTURE &denominator) {
    using CR = Process_Amp::CTH_RES;
    CR n0 = numerator[0];
    CR n1 = numerator[1];
    CR n2 = numerator[2];

    using PS = Process_Amp::PROPAGATOR_STRUCTURE;
    PS p0 = denominator.first;
    PS p1 = denominator.second;

    int id_p0 = p0.first;
    CR d00 = p0.second[0];
    CR d01 = p0.second[1];

    int id_p1 = p1.first;
    CR d10 = p1.second[0];
    CR d11 = p1.second[1];

    REAL nn0 = n0.second;
    REAL nn1 = n1.second;
    REAL nn2 = n2.second;
    REAL nd00 = d00.second;
    REAL nd01 = d01.second;
    REAL nd10 = d10.second;
    REAL nd11 = d11.second;

    // * The amplitude has following structure
    // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10+d11*cth)
    // * Then according to different situation for {d00,d01} and {d10,d11}, we will formally integral over cth
    if (id_p0 == id_p1) {
        // * In this case d00==d10, and d01 == d11
        // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)^2
        if (d01.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)^2
            return INT_RES_11_1(nn0, nn1, nn2, nd00, nd01);
        }

        if (!d01.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00)^2
            return INT_RES_00(nn0, nn1, nn2, nd00, nd00);
        }
    }

    if (id_p0 != id_p1) {
        // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10 + d11*cth)
        if (d01.first && d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10 + d11*cth)
            return INT_RES_11_0(nn0, nn1, nn2, nd00, nd01, nd10, nd11);
        }

        if (d01.first && !d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10)
            return INT_RES_01(nn0, nn1, nn2, nd00, nd01, nd10);
        }

        if (!d01.first && d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00)/(d10 + d11*cth)
            return INT_RES_01(nn0, nn1, nn2, nd10, nd11, nd00);
        }

        if (!d01.first && !d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00)/(d10)
            return INT_RES_00(nn0, nn1, nn2, nd00, nd10);
        }
    }

    return 0;
}
