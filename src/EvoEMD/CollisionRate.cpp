#include "EvoEMD/CollisionRate.h"

#include <cmath>

#include "EvoEMD/Constants.h"
#include "EvoEMD/PhaseSpace.h"
#include "EvoEMD/spdlog_wrapper.h"
#include "cuba.h"
#include "gsl/gsl_sf_bessel.h"

namespace EvoEMD {

// REAL Decay12_Rate::Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) {
//     // * For  1 -> 2 3;
//     REAL m1 = amp->INITIAL[0]->Get_Mass();
//     REAL m2 = amp->FINAL[0]->Get_Mass();
//     REAL m3 = amp->FINAL[1]->Get_Mass();
//     REAL m12 = m1 * m1;
//     REAL m22 = m2 * m2;
//     REAL m32 = m3 * m3;

//     // * Final state momentum when 1 is at rest;
//     REAL pf = sqrt(Kallen_Lam(m12, m22, m32)) / 2.0 / m1;

//     // * Two body phase space volumn
//     const Process_Amp &res = amp->Get_Amp(sqrt_shat);
//     REAL PS2 = 1.0 / (2.0 * m1) * pf / (16.0 * M_PI * M_PI * m1);
//     REAL RES = 2.0                      // * From integral over cos(theta)
//                * 2.0 * M_PI             // * From integral over phi
//                * (res[0].Numerator[0])  // * The amplitude
//                * PS2;                   // * The two body phase space
//     return RES;
// }

REAL Decay12_Rate::Get_Collision_Rate(REAL T) {
    REAL AMP_WITH_SOLID_ANGLE = amp->Get_Amp(0);
    REAL m1 = amp->INITIAL[0]->Get_Mass();
    REAL m2 = amp->FINAL[0]->Get_Mass();
    REAL m3 = amp->FINAL[1]->Get_Mass();
    REAL m12 = m1 * m1;
    REAL m22 = m2 * m2;
    REAL m32 = m3 * m3;

    REAL pf = sqrt(Kallen_Lam(m12, m22, m32)) / 2.0 / m1;
    REAL PS2 = 1.0 / (2.0 * m1) * pf / (16.0 * M_PI * M_PI * m1);

    REAL z = m1 / T;
    REAL k1z = z > BESSEL_Z_MAX ? 0 : gsl_sf_bessel_K1(z);
    REAL RES = m1 * m1 * T / 2.0 / M_PI / M_PI * k1z * PS2 * AMP_WITH_SOLID_ANGLE;
    return RES;
}

// REAL Scatter22_Rate::Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) {
//     const Process_Amp &amp_res = amp->Get_Amp(sqrt_shat);
//     REAL amp_total = 0;
//     for (unsigned i_diag = 0; i_diag < amp_res.size(); i_diag++) {
//         amp_total += Get_Amp_Integrate_over_Phase_Space_Single_Channel(amp_res[i_diag]);
//     }
//     return amp_total;  // * Extra 2pi from integrate over phi;
// }

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

// REAL Scatter22_Rate::Get_Amp_Integrate_over_Phase_Space_Single_Channel(const Process_Amp_Single_Diagram &amp_single)
// {
//     using NT = Process_Amp_Single_Diagram::Numerator_Type;
//     using DT = Process_Amp_Single_Diagram::Denominator_Type;
//     REAL n0 = amp_single.Numerator[0];
//     REAL n1 = amp_single.Numerator[1];
//     REAL n2 = amp_single.Numerator[2];

//     using PT = Process_Amp_Single_Diagram::Propagator_Type;
//     using RT = Process_Amp_Single_Diagram::Result_Type;
//     const DT &denominator = amp_single.Denominator;
//     const PT &p0 = denominator.first;
//     const PT &p1 = denominator.second;

//     int id_p0 = p0.first;
//     const RT &crf_p0 = p0.second;
//     int ncth_p0 = crf_p0.size();

//     int id_p1 = p1.first;
//     const RT &crf_p1 = p1.second;
//     int ncth_p1 = crf_p1.size();

//     // * The amplitude has following structure
//     // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10+d11*cth)
//     // * Then according to different situation for {d00,d01} and {d10,d11}, we will formally integral over cth
//     if (id_p0 == id_p1) {
//         // * In this case d00==d10, and d01 == d11
//         // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)^2
//         if (ncth_p0 == 2) {
//             // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)^2
//             // * 2pi from integral over phi, INT_RES_XXXX from integral over cth
//             REAL d00 = crf_p0[0];
//             REAL d01 = crf_p0[1];
//             return 2.0 * M_PI * INT_RES_11_1(n0, n1, n2, d00, d01);
//         }

//         if (ncth_p0 == 1) {
//             // * (n0 + n1*cth + n2*cth^2)/(d00)^2
//             REAL d00 = crf_p0[0];
//             return 2.0 * M_PI * INT_RES_00(n0, n1, n2, d00, d00);
//         }
//     }

//     if (id_p0 != id_p1) {
//         // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10 + d11*cth)
//         if (ncth_p0 == 2 && ncth_p1 == 2) {
//             // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10 + d11*cth)
//             REAL d00 = crf_p0[0];
//             REAL d01 = crf_p0[1];
//             REAL d10 = crf_p1[0];
//             REAL d11 = crf_p1[1];
//             return 2.0 * M_PI * INT_RES_11_0(n0, n1, n2, d00, d01, d10, d11);
//         }

//         if (ncth_p0 == 2 && ncth_p1 == 1) {
//             // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10)
//             REAL d00 = crf_p0[0];
//             REAL d01 = crf_p0[1];
//             REAL d10 = crf_p1[0];
//             return 2.0 * M_PI * INT_RES_01(n0, n1, n2, d00, d01, d10);
//         }

//         if (ncth_p0 == 1 && ncth_p1 == 2) {
//             // * (n0 + n1*cth + n2*cth^2)/(d00)/(d10 + d11*cth)
//             REAL d00 = crf_p0[0];
//             REAL d10 = crf_p1[0];
//             REAL d11 = crf_p1[1];
//             return 2.0 * M_PI * INT_RES_01(n0, n1, n2, d10, d11, d00);
//         }

//         if (ncth_p0 == 1 && ncth_p1 == 1) {
//             // * (n0 + n1*cth + n2*cth^2)/(d00)/(d10)
//             REAL d00 = crf_p0[0];
//             REAL d10 = crf_p1[0];
//             return 2.0 * M_PI * INT_RES_00(n0, n1, n2, d00, d10);
//         }
//     }

//     return 0;
// }

struct INT_PARAM {
    Scatter22_Rate *ptr;
    REAL Temperature;
    REAL sqrt_s_min;
    REAL sqrt_s_max;
    REAL m1;
    REAL m2;
    REAL m3;
    REAL m4;
};
int Scatter22_Rate_Integrand(const int *ndim, const REAL x[], const int *ncomp, REAL ff[], void *params) {
    INT_PARAM *par = (INT_PARAM *)params;
    REAL T = par->Temperature;

    REAL sqrt_s_min = par->sqrt_s_min;
    REAL sqrt_s_max = par->sqrt_s_max;
    REAL sqrt_s = sqrt_s_min * exp(x[0] * log(sqrt_s_max / sqrt_s_min));
    REAL Jac_from_sqrt_s = log(sqrt_s_max / sqrt_s_min) * sqrt_s;
    REAL s = sqrt_s * sqrt_s;
    REAL m1 = par->m1;
    REAL m2 = par->m2;
    REAL m3 = par->m3;
    REAL m4 = par->m4;

    REAL sqrt_Kallen_init = sqrt(Kallen_Lam(1.0, m1 * m1 / s, m2 * m2 / s));
    REAL sqrt_Kallen_final = sqrt(Kallen_Lam(1.0, m3 * m3 / s, m4 * m4 / s));

    REAL k1z = sqrt_s / T > BESSEL_Z_MAX ? 0 : gsl_sf_bessel_K1(sqrt_s / T);
    ff[0] =
        Jac_from_sqrt_s * s * k1z * sqrt_Kallen_init * sqrt_Kallen_final * par->ptr->Get_Amp_With_Solid_Angle(sqrt_s);
    return 0;
}

// CUBA PARAMETERS:
// COMMON VARIABLES
#define MINEVAL 5000
#define MAXEVAL 100000
#define EPSREL 1E-3
#define EPSABS 0  // 1E-7 precision is about 1 ab, which is way small.
#define NVEC 1
#define FLAGS 0
#define SEED 0

// VEGAS
#define NSTART 2000
#define NINCREASE 500
#define NBATCH 2000
#define GRIDNO 0

REAL Scatter22_Rate::Get_Collision_Rate(REAL T) {
    REAL m1 = amp->INITIAL[0]->Get_Mass();
    REAL m2 = amp->INITIAL[1]->Get_Mass();
    REAL m3 = amp->FINAL[0]->Get_Mass();
    REAL m4 = amp->FINAL[1]->Get_Mass();

    INT_PARAM par;
    par.ptr = this;
    par.Temperature = T;
    par.sqrt_s_min = std::max(m1 + m2, m3 + m4);
    par.sqrt_s_max = std::max(100 * par.sqrt_s_min, 100 * T);
    par.m1 = m1;
    par.m2 = m2;
    par.m3 = m3;
    par.m4 = m4;

    int fail, neval;
    const int NDIM = 1;
    const int NCOMP = 1;
    REAL RES[NCOMP];
    REAL ERR[NCOMP];
    REAL PROB[NCOMP];

    Vegas(NDIM, NCOMP, Scatter22_Rate_Integrand, &par, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART,
          NINCREASE, NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);

    RES[0] *= T / 16.0 / pow(2.0 * M_PI, 6);
    return RES[0];
}

}  // namespace EvoEMD
