#include "EvoEMD/CollisionRate.h"

#include <cmath>

#include "EvoEMD/Constants.h"
#include "EvoEMD/PhaseSpace.h"
#include "EvoEMD/spdlog_wrapper.h"
#include "cuba.h"
#include "gsl/gsl_sf_bessel.h"

namespace EvoEMD {

REAL Decay_Rate::Get_Collision_Rate(REAL T) {
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

struct INT_PARAM {
    Scatter_Rate *ptr;
    REAL Temperature;
    REAL sqrt_s_min;
    REAL sqrt_s_max;
    REAL m1;
    REAL m2;
    REAL m3;
    REAL m4;
};
int Scatter_Rate_Integrand(const int *ndim, const REAL x[], const int *ncomp, REAL ff[], void *params) {
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

REAL Scatter_Rate::Get_Collision_Rate(REAL T) {
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

    Vegas(NDIM, NCOMP, Scatter_Rate_Integrand, &par, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART,
          NINCREASE, NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);

    RES[0] *= T / 16.0 / pow(2.0 * M_PI, 6);
    return RES[0];
}

}  // namespace EvoEMD
