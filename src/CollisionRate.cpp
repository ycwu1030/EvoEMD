#include "CollisionRate.h"

#include <cmath>

#include "PhaseSpace.h"
#include "cuba.h"
#include "gsl/gsl_sf_bessel.h"

const Process_Amp::NUMERATOR_STRUCTURE &Process_Amp::Get_Numerator(int diagram_id) const {
    return amps_numerator.at(diagram_id);
    // return res;
}

const Process_Amp::DENOMINATOR_STRUCTURE &Process_Amp::Get_Denominator(int diagram_id) const {
    return amps_denominator.at(diagram_id);
}

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
    REAL k1z = gsl_sf_bessel_K1(z);
    REAL RES = m1 * m1 * T / 2.0 / M_PI / M_PI * k1z;
    return RES;
}

REAL Scatter22_Rate::Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) {
    using NS = Process_Amp::NUMERATOR_STRUCTURE;
    using DS = Process_Amp::DENOMINATOR_STRUCTURE;
    const Process_Amp &amp_res = amp->Get_Amp(sqrt_shat);
    REAL amp_total = 0;
    for (unsigned i_diag = 0; i_diag < amp_res.n_diag; i_diag++) {
        const NS &num = amp_res.Get_Numerator(i_diag);
        const DS &den = amp_res.Get_Denominator(i_diag);
        amp_total += Get_Amp_Integrate_over_Phase_Space_Single_Channel(num, den);
    }
    return amp_total;  // * Extra 2pi from integrate over phi;
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
            // * 2pi from integral over phi, INT_RES_XXXX from integral over cth
            return 2.0 * M_PI * INT_RES_11_1(nn0, nn1, nn2, nd00, nd01);
        }

        if (!d01.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00)^2
            return 2.0 * M_PI * INT_RES_00(nn0, nn1, nn2, nd00, nd00);
        }
    }

    if (id_p0 != id_p1) {
        // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10 + d11*cth)
        if (d01.first && d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10 + d11*cth)
            return 2.0 * M_PI * INT_RES_11_0(nn0, nn1, nn2, nd00, nd01, nd10, nd11);
        }

        if (d01.first && !d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00 + d01*cth)/(d10)
            return 2.0 * M_PI * INT_RES_01(nn0, nn1, nn2, nd00, nd01, nd10);
        }

        if (!d01.first && d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00)/(d10 + d11*cth)
            return 2.0 * M_PI * INT_RES_01(nn0, nn1, nn2, nd10, nd11, nd00);
        }

        if (!d01.first && !d11.first) {
            // * (n0 + n1*cth + n2*cth^2)/(d00)/(d10)
            return 2.0 * M_PI * INT_RES_00(nn0, nn1, nn2, nd00, nd10);
        }
    }

    return 0;
}

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

    REAL k1z = gsl_sf_bessel_K1(sqrt_s / T);
    ff[0] = Jac_from_sqrt_s * s * k1z * sqrt_Kallen_init * sqrt_Kallen_final *
            par->ptr->Get_Amp_Integrate_over_Phase_Space(sqrt_s);
    return 0;
}

// CUBA PARAMETERS:
// COMMON VARIABLES
#define MINEVAL 10000
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

Process::Process(Amplitude *amp) {
    this->amp = amp;
    INIT = amp->INITIAL;
    FINAL = amp->FINAL;
    if (amp->N_INITIAL == 1 && amp->N_FINAL == 2) {
        CR_Calculator = new Decay12_Rate(amp);
    }
    if (amp->N_INITIAL == 2 && amp->N_FINAL == 2) {
        CR_Calculator = new Scatter22_Rate(amp);
    }
    for (int i = 0; i < INIT.size(); i++) {
        INIT[i]->Register_Process(this);
    }
    for (int i = 0; i < FINAL.size(); i++) {
        FINAL[i]->Register_Process(this);
    }
}

Process::~Process() { delete CR_Calculator; }

REAL Process::Get_Collision_Rate(REAL T) { return CR_Calculator->Get_Collision_Rate(T); }

REAL Process::Get_Yield_Coeff(REAL T) { return amp->Get_Coeff(T); }
