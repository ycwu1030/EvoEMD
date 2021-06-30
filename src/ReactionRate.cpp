#include <algorithm>
#include <cmath>

#include "LeptogenesisRate.h"
#include "PhaseSpace.h"
#include "cuba.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_sf_bessel.h"
#include "spdlog_wrapper.h"

using namespace std;

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

REAL LeptogenesisRate::Calc_NChiS_Gamma(REAL Temp, int i) {
    SPDLOG_INFO_FILE("Calculate gammas for N{} -> Chi S at T = {:+9.8e}", i + 1, Temp);
    REAL MN = i == 0 ? MNR1 : MNR2;
    if (MN < MCHI + MS) {
        SPDLOG_INFO_FILE("Not Sufficient Phase Space, Return 0.");
        return 0;
    }
    REAL A = pow(MN / MNR1, 2);
    REAL B = pow(MCHI / MNR1, 2);
    REAL C = pow(MS / MNR1, 2);
    REAL z = MNR1 / Temp;
    REAL SQLam = sqrt(Kallen_Lam(1, B / A, C / A));
    REAL gamma = LamX * LamX / 16.0 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) /
                 z * (1 + B / A - C / A) * SQLam;
    SPDLOG_INFO_FILE("gamma(N{} -> Chi S) = {:+9.8e}", i + 1, gamma);
    return gamma;
}

REAL LeptogenesisRate::Calc_ChiSN_Gamma(REAL Temp, int i) {
    SPDLOG_INFO_FILE("Calculate gammas for Chi -> S N{} at T = {:+9.8e}", i + 1, Temp);
    REAL MN = i == 0 ? MNR1 : MNR2;
    if (MCHI < MN + MS) {
        SPDLOG_INFO_FILE("Insufficient Phase Space, Return 0.");
        return 0;
    }
    REAL A = pow(MCHI / MNR1, 2);
    REAL B = pow(MN / MNR1, 2);
    REAL C = pow(MS / MNR1, 2);
    REAL z = MNR1 / Temp;
    REAL SQLam = sqrt(Kallen_Lam(1, B / A, C / A));
    REAL gamma = LamX * LamX / 8.0 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) / z *
                 (1 + B / A - C / A) * SQLam;
    SPDLOG_INFO_FILE("gamma(Chi -> S N{}) = {:+9.8e}", i + 1, gamma);
    return gamma;
}

REAL LeptogenesisRate::Calc_SChiN_Gamma(REAL Temp, int i) {
    SPDLOG_INFO_FILE("Calculate gammas for S -> Chi N{} at T = {:+9.8e}", i + 1, Temp);
    REAL MN = i == 0 ? MNR1 : MNR2;
    if (MS < MN + MCHI) {
        SPDLOG_INFO_FILE("Insufficient Phase Space, Return 0.");
        return 0;
    }
    REAL A = pow(MS / MNR1, 2);
    REAL B = pow(MN / MNR1, 2);
    REAL C = pow(MCHI / MNR1, 2);
    REAL z = MNR1 / Temp;
    REAL SQLam = sqrt(Kallen_Lam(1, B / A, C / A));
    REAL gamma = LamX * LamX / 16.0 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) /
                 z * (1 - B / A - C / A) * SQLam;
    SPDLOG_INFO_FILE("gamma(S -> Chi N{}) = {:+9.8e}.", i + 1, gamma);
    return gamma;
}

REAL LeptogenesisRate::Calc_NLPhi_Gamma(REAL Temp, int i) {
    SPDLOG_INFO_FILE("Calculate gammas for N{} -> L Phi at T = {:+9.8e}", i + 1, Temp);
    REAL z = MNR1 / Temp;
    REAL MN = i == 0 ? MNR1 : MNR2;
    REAL A = pow(MN / MNR1,2);
    REAL coup = real(Nu_Param.Get_YdagYij(i, i));
    REAL gamma = coup * coup / 8 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) / z;
    SPDLOG_INFO_FILE("gamma(N{} -> L Phi) = {:+9.8e}.", i + 1, gamma);
    return gamma;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen(REAL s, int processid, int i, int j) {
    switch (processid) {
        case 0:
            return SqAmp_dOmega_with_Kallen_LPhiChiS(s);
        case 1:
            return SqAmp_dOmega_with_Kallen_NNChiChi(s, i, j);
        case 2:
            return SqAmp_dOmega_with_Kallen_NNSS(s, i, j);
        default:
            return SqAmp_dOmega_with_Kallen_LPhiChiS(s);
    }
}

struct INTE_PARAM {
    LeptogenesisRate *ptr;
    int processid;
    int ni;
    int nj;
    REAL Temp;
    REAL smin;
    REAL smax;
};

int gamma_Integrand(const int *ndim, const REAL x[], const int *ncomp, REAL ff[], void *_params) {
    INTE_PARAM *params = (INTE_PARAM *)_params;
    int procid = params->processid;
    int ni = params->ni;
    int nj = params->nj;
    REAL Temp = params->Temp;
    REAL smin = params->smin;
    REAL smax = params->smax;
    REAL s = smin + (smax - smin) * x[0];
    REAL JacS = (smax - smin); 
    // REAL s = smin*exp(x[0]*log(smax/smin));
    // JacS = smin*log(smax/smin)*exp(x[0]*log(smax/smin));


    ff[0] = JacS * params->ptr->SqAmp_dOmega_with_Kallen(s, procid, ni, nj) * sqrt(s) *
            gsl_sf_bessel_K1(sqrt(s) / Temp) * Temp / 32.0 / pow(2 * M_PI, 6);
    return 0;
}

REAL LeptogenesisRate::Calc_LPhiChiS_Gamma(REAL Temp) {
    SPDLOG_INFO_FILE("Calculate gammas for L Phi <-> Chi S at T = {:+9.8e}.", Temp);
    INTE_PARAM param = {this, 0, -1, -1, Temp, pow(MCHI + MS, 2), pow(100 * Temp, 2)};
    int fail, neval;
    const int NDIM = 1;
    const int NCOMP = 1;
    REAL RES[NCOMP];
    REAL ERR[NCOMP];
    REAL PROB[NCOMP];

    Vegas(NDIM, NCOMP, gamma_Integrand, &param, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART, NINCREASE,
          NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);
    SPDLOG_INFO_FILE("gamma(L Phi <-> Chi S) = {:+9.8e}.", RES[0]);
    return RES[0];
}

REAL LeptogenesisRate::Calc_NNChiChi_Gamma(REAL Temp, int i, int j) {
    SPDLOG_INFO_FILE("Calculate gammas for N{} N{} <-> Chi Chi at T = {:+9.8e}.", i + 1, j + 1, Temp);
    REAL MNI = i == 0 ? MNR1 : MNR2;
    REAL MNJ = j == 0 ? MNR1 : MNR2;
    REAL mmin = max(MNI + MNJ, 2 * MCHI);
    INTE_PARAM param = {this, 1, i, j, Temp, pow(mmin, 2), pow(100 * Temp, 2)};
    int fail, neval;
    const int NDIM = 1;
    const int NCOMP = 1;
    REAL RES[NCOMP];
    REAL ERR[NCOMP];
    REAL PROB[NCOMP];

    Vegas(NDIM, NCOMP, gamma_Integrand, &param, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART, NINCREASE,
          NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);

    SPDLOG_INFO_FILE("gamma(N{} N{} <-> Chi Chi) = {:+9.8e}", i + 1, j + 1, RES[0]);
    return RES[0];
}

REAL LeptogenesisRate::Calc_NNSS_Gamma(REAL Temp, int i, int j) {
    SPDLOG_INFO_FILE("Calculate gammas for N{} N{} <-> S S at T = {:+9.8e}.", i + 1, j + 1, Temp);
    REAL MNI = i == 0 ? MNR1 : MNR2;
    REAL MNJ = j == 0 ? MNR1 : MNR2;
    REAL mmin = max(MNI + MNJ, 2 * MCHI);
    INTE_PARAM param = {this, 2, i, j, Temp, pow(mmin, 2), pow(100 * Temp, 2)};
    int fail, neval;
    const int NDIM = 1;
    const int NCOMP = 1;
    REAL RES[NCOMP];
    REAL ERR[NCOMP];
    REAL PROB[NCOMP];

    Vegas(NDIM, NCOMP, gamma_Integrand, &param, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART, NINCREASE,
          NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);

    SPDLOG_INFO_FILE("gamma(N{} N{} <-> S S) = {:+9.8e}.", i + 1, j + 1, Temp);
    return RES[0];
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen_LPhiChiS(REAL s) {
    REAL M1 = 0;
    REAL M2 = 0;
    REAL M3 = MCHI;
    REAL M4 = MS;
    if (s < pow(M3 + M4, 2) || s < pow(M1 + M2, 2)) {
        return 0;
    }
    REAL SQLamFactor = sqrt(Kallen_Lam(1.0, M1 * M1 / s, M2 * M2 / s) * Kallen_Lam(1.0, M3 * M3 / s, M4 * M4 / s));
    REAL prefix = 2 * M_PI * LamX * LamX * (s + MCHI * MCHI - MS * MS);
    REAL h1a2coup = 0;
    REAL h2a2coup = 0;
    REAL Rh1aha2coup = 0;
    REAL Ih1aha2coup = 0;
    for (int lep_index = 0; lep_index < 3; lep_index++) {
        auto h1a = Nu_Param.Get_Yij(0, lep_index);
        auto h2a = Nu_Param.Get_Yij(1, lep_index);
        h1a2coup += norm(h1a);  // norm is the square of its magnitude
        h2a2coup += norm(h2a);
        Rh1aha2coup += real(h1a * conj(h2a));
        Ih1aha2coup += imag(h1a * conj(h2a));
    }
    REAL N1Diag = (s + MNR1 * MNR1) / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2));
    REAL N2Diag = (s + MNR2 * MNR2) / (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterPref = 2 * (s + MNR1 * MNR2) / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2)) /
                     (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterR = InterPref * ((s - MNR1 * MNR1) * (s - MNR2 * MNR2) + MNR1 * MNR2 * GammaN1 * GammaN2);
    REAL InterI = InterPref * ((s - MNR1 * MNR1) * MNR2 * GammaN2 - (s - MNR2 * MNR2) * MNR1 * GammaN1);
    REAL RES =
        SQLamFactor * prefix * (h1a2coup * N1Diag + h2a2coup * N2Diag + Rh1aha2coup * InterR + Ih1aha2coup * InterI);
    return RES;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen_LPhiLbarPhiBar(REAL s) {
    if (s < 0) {
        return 0;
    }
    REAL SQLamFactor = 1.;
    REAL prefix = 2 * M_PI * s*s;
    REAL h1a2coup = 0;
    REAL h2a2coup = 0;
    REAL h1b2coup = 0;
    REAL h2b2coup = 0;
    REAL Rh1aha2h1bhb2coup = 0;
    REAL Ih1aha2h1bhb2coup = 0;
    for (int lep_index_a = 0; lep_index_a < 3; lep_index_a++) {
        for (int lep_index_b = 0; lep_index_b < 3; lep_index_b++) {
            auto h1a = Nu_Param.Get_Yij(0, lep_index_a);
            auto h2a = Nu_Param.Get_Yij(1, lep_index_a);
            auto h1b = Nu_Param.Get_Yij(0, lep_index_b);
            auto h2b = Nu_Param.Get_Yij(1, lep_index_b);
            h1a2coup += norm(h1a);  // norm is the square of its magnitude
            h1b2coup += norm(h1b);  
            h2a2coup += norm(h2a);
            h2b2coup += norm(h2b);
            Rh1aha2h1bhb2coup += real(h1a * conj(h2a) * h1b * conj(h2b));
            Ih1aha2h1bhb2coup += imag(h1a * conj(h2a) * h1b * conj(h2b));
        }
    }
    REAL N1Diag = 1. / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2));
    REAL N2Diag = 1. / (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterPref = 2. / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2)) / (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterR = InterPref * ((s - MNR1 * MNR1) * (s - MNR2 * MNR2) + MNR1 * MNR2 * GammaN1 * GammaN2);
    REAL InterI = InterPref * ((s - MNR1 * MNR1) * MNR2 * GammaN2 - (s - MNR2 * MNR2) * MNR1 * GammaN1);
    REAL RES =
        SQLamFactor * prefix * (h1a2coup * h1b2coup * N1Diag + h2a2coup * h2b2coup * N2Diag + Rh1aha2h1bhb2coup * InterR + Ih1aha2h1bhb2coup * InterI);
    return RES;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen_NNChiChi(REAL s, int i, int j) {
    REAL MNI = i == 0 ? MNR1 : MNR2;
    REAL MNJ = j == 0 ? MNR1 : MNR2;
    if (s < pow(MNI + MNJ, 2) || s < pow(2 * MCHI, 2)) {
        return 0;
    }
    REAL SQLamFactor =
        sqrt(Kallen_Lam(1.0, MNI * MNI / s, MNJ * MNJ / s) * Kallen_Lam(1.0, MCHI * MCHI / s, MCHI * MCHI / s));
    REAL prefix = 2 * M_PI * pow(LamX, 4);
    REAL EI = Ei(sqrt(s), MNI, MNJ);
    REAL EJ = Ei(sqrt(s), MNJ, MNI);
    REAL Pc = Pi(sqrt(s), MNI, MNJ);
    REAL E1 = Ei(sqrt(s), MCHI, MCHI);
    REAL E2 = Ei(sqrt(s), MCHI, MCHI);
    REAL Pf = Pi(sqrt(s), MCHI, MCHI);
    // * T-channel
    REAL a = E1 * EI;
    REAL b = E2 * EJ;
    REAL c = Pc * Pf;
    REAL d = 2 * E1 * EI + MS * MS - MNI * MNI - MCHI * MCHI;
    // * Int[(a - c x)(b - c x)/(d - 2 c x)/(d - 2 c x),{x, -1, 1}]
    REAL RES_T = 4 * (1.0 / 2.0 - (d - a - b) / 4 / c * log((d + 2 * c) / (d - 2 * c)) +
                      (d - 2 * a) * (d - 2 * b) / 2 / (d - 2 * c) / (d + 2 * c));
    // * U-channel
    a = E2 * EI;
    b = E1 * EJ;
    c = Pc * Pf;
    d = 2 * E1 * EJ + MS * MS - MNJ * MNJ - MCHI * MCHI;
    // * Int[(a + c x)(b + c x)/(d + 2 c x)/(d + 2 c x),{x, -1, 1}]
    REAL RES_U = 4 * (1.0 / 2.0 - (d - a - b) / 4 / c * log((d + 2 * c) / (d - 2 * c)) +
                      (d - 2 * a) * (d - 2 * b) / 2 / (d - 2 * c) / (d + 2 * c));
    // * Interference between t and u channel
    a = 2 * E1 * EI + MS * MS - MNI * MNI - MCHI * MCHI;
    b = 2 * E1 * EJ + MS * MS - MNJ * MNJ - MCHI * MCHI;
    c = Pc * Pf;
    // * Int[1/(a - 2 c x)/(b + 2 c x), {x, -1 ,1}]
    REAL RES_INTER = 2 * MNI * MNJ * (2 * MCHI * MCHI - s) *
                     (log((a + 2 * c) / (a - 2 * c)) + log((b + 2 * c) / (b - 2 * c))) / 2 / (a + b) / c;
    REAL RES = prefix * SQLamFactor * (RES_T + RES_U + RES_INTER);
    return RES;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen_NNSS(REAL s, int i, int j) {
    REAL MNI = i == 0 ? MNR1 : MNR2;
    REAL MNJ = j == 0 ? MNR1 : MNR2;
    if (s < pow(MNI + MNJ, 2) || s < pow(2 * MS, 2)) {
        return 0;
    }
    REAL SQLamFactor = sqrt(Kallen_Lam(1.0, MNI * MNI / s, MNJ * MNJ / s) * Kallen_Lam(1.0, MS * MS / s, MS * MS / s));
    REAL prefix = 2 * M_PI * pow(LamX, 4);
    REAL EI = Ei(sqrt(s), MNI, MNJ);
    REAL EJ = Ei(sqrt(s), MNJ, MNI);
    REAL Pc = Pi(sqrt(s), MNI, MNJ);
    REAL E1 = Ei(sqrt(s), MS, MS);
    REAL E2 = Ei(sqrt(s), MS, MS);
    REAL Pf = Pi(sqrt(s), MS, MS);

    REAL c = 2 * Pc * Pf;
    REAL t_neg_a = 2 * E1 * EI - MNI * MNI - MS * MS;
    REAL u_neg_a = 2 * E1 * EJ - MNJ * MNJ - MS * MS;
    // * T-Channel
    // TODO: write the expression into the note
    REAL RES_T =
        (MS * MS - MNI * MNI) * (MS * MS - MNJ * MNJ) * 2 / (pow(t_neg_a + MCHI * MCHI, 2) - pow(c, 2)) -
        (MNI * MNI + MNJ * MNJ) * (log((t_neg_a + MCHI * MCHI + c) / (t_neg_a + MCHI * MCHI - c)) / c -
                                   2 * MCHI * MCHI / (pow(t_neg_a + MCHI * MCHI, 2) - pow(c, 2))) +
        (2 -
         (t_neg_a + u_neg_a + 2 * MCHI * MCHI) * log((t_neg_a + MCHI * MCHI + c) / (t_neg_a + MCHI * MCHI - c)) / c +
         2 * (t_neg_a + u_neg_a + MCHI * MCHI) * (MCHI * MCHI) / (pow(t_neg_a + MCHI * MCHI, 2) - pow(c, 2)));
    // * U-Channel
    REAL RES_U =
        (MS * MS - MNI * MNI) * (MS * MS - MNJ * MNJ) * 2 / (pow(u_neg_a + MCHI + MCHI, 2) - pow(c, 2)) -
        (MNI * MNI + MNJ * MNJ) * (log((u_neg_a + MCHI + MCHI + c) / (u_neg_a + MCHI + MCHI - c)) / c -
                                   2 * MCHI * MCHI / (pow(u_neg_a + MCHI + MCHI, 2) - pow(c, 2))) +
        (2 -
         (u_neg_a + t_neg_a + 2 * MCHI * MCHI) * log((u_neg_a + MCHI + MCHI + c) / (u_neg_a + MCHI + MCHI - c)) / c +
         2 * (u_neg_a + t_neg_a + MCHI * MCHI) * (MCHI * MCHI) / (pow(u_neg_a + MCHI + MCHI, 2) - pow(c, 2)));
    // * Interference between t and u channel
    REAL RES_INTER = 2 * MNI * MNI * MNJ * MNJ * (MNI * MNI + MNJ * MNJ - 2 * MS * MS) / (u_neg_a + t_neg_a) / c *
                     (log((u_neg_a + MCHI + MCHI + c) / (u_neg_a + MCHI + MCHI - c)) +
                      log((t_neg_a + MCHI * MCHI + c) / (t_neg_a + MCHI * MCHI - c)));
    return prefix * SQLamFactor * (RES_T + RES_U + RES_INTER);
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen_NNChiChi_v2(REAL s, int i, int j) {
    REAL MNI = i == 0 ? MNR1 : MNR2;
    REAL MNJ = j == 0 ? MNR1 : MNR2;
    if (s < pow(MNI + MNJ, 2) || s < pow(2 * MCHI, 2)) {
        return 0;
    }
    REAL SQLamFactor =
        sqrt(Kallen_Lam(1.0, MNI * MNI / s, MNJ * MNJ / s) * Kallen_Lam(1.0, MCHI * MCHI / s, MCHI * MCHI / s));
    REAL prefix = 2 * M_PI * pow(LamX, 4);
    
    REAL E1 = Ei(sqrt(s), MNI, MNJ);
    REAL E2 = Ei(sqrt(s), MNJ, MNI);
    REAL E3 = Ei(sqrt(s), MCHI, MCHI);
    REAL E4 = Ei(sqrt(s), MCHI, MCHI);
    REAL P12 = Pi(sqrt(s), MNI, MNJ);
    REAL P34 = Pi(sqrt(s), MCHI, MCHI);
    
    REAL a = 4 * E1 * E3;
    REAL b = 4 * E2 * E4;
    REAL c = 4 * P12 * P34;
    REAL d = 2 * MS * MS + 4 * E1 * E3 - 2 * MNI * MNI - 2 * MCHI * MCHI;
    // Conditional expression for the integration over cosine
    if (d < c) {
        return 0;
    }

    // * T-channel
    // * Int[(a - c x)(b - c x)/(d - c x)/(d - c x),{x, -1, 1}]
    REAL RES_T = (a + b - 2 * d)/c * log((d + c) / (d - c)) -
                      2 * (c*c + d * (b - 2 * d) + a * (d - b)) / (d*d - c*c);
    // * U-channel
    // * Int[(a + c x)(b + c x)/(d + c x)/(d + c x),{x, -1, 1}]
    REAL RES_U = RES_T;
    // * Interference between t and u channel
    // * propto Int[1/(d^2 - c^2 x), {x, -1 ,1}]
    REAL RES_INTER = 8 * MNI * MNJ * (s - 2 * MCHI * MCHI) / c / d * log((d + c) / (d - c));
    REAL RES = prefix * SQLamFactor * (RES_T + RES_U + RES_INTER);
    return RES;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen_NNSS_v2(REAL s, int i, int j) {
    REAL MNI = i == 0 ? MNR1 : MNR2;
    REAL MNJ = j == 0 ? MNR1 : MNR2;
    if (s < pow(MNI + MNJ, 2) || s < pow(2 * MS, 2)) {
        return 0;
    }
    REAL SQLamFactor =
        sqrt(Kallen_Lam(1.0, MNI * MNI / s, MNJ * MNJ / s) * Kallen_Lam(1.0, MS * MS / s, MS * MS / s));
    REAL prefix = 2 * M_PI * pow(LamX, 4);
    
    REAL E1 = Ei(sqrt(s), MNI, MNJ);
    REAL E2 = Ei(sqrt(s), MNJ, MNI);
    REAL E3 = Ei(sqrt(s), MS, MS);
    REAL E4 = Ei(sqrt(s), MS, MS);
    REAL P12 = Pi(sqrt(s), MNI, MNJ);
    REAL P34 = Pi(sqrt(s), MS, MS);
    
    REAL a = 4 * E1 * E3;
    REAL b = 4 * E2 * E4;
    REAL c = 4 * P12 * P34;
    REAL d = 2 * MCHI*MCHI + 4 * E1*E3 - 2 * MNI*MNI - 2 * MS*MS;
    REAL g = 2 * ( (s - MNJ*MNJ) * (s - 4*MS*MS + MNJ*MNJ) - pow(MNI,4) - MNI*MNI * (6 * MNJ*MNJ - 4 * MS*MS )  );
    // Conditional expression for the integration over cosine
    if (d < c) {
        return 0;
    }

    // * T-channel
    // * Int[(g - 4(m1^2 + m2^2) c x - 2c^2 x^2)/2/(d - c x)/(d - c x),{x, -1, 1}]
    REAL RES_T = 2 * (d + MNI*MNI + MNJ*MNJ)/c * log((d + c) / (d - c)) +
                      (2 * c*c + g - 4 * d * (d + MNI*MNI + MNJ*MNJ)) / (d*d - c*c);
    // * U-channel
    // * Int[(g + 4(m1^2 + m2^2) c x - 2c^2 x^2)/2/(d + c x)/(d + c x),{x, -1, 1}]
    REAL RES_U = RES_T;
    // * Interference between t and u channel
    // * propto Int[1/(d^2 - c^2 x), {x, -1 ,1}]
    REAL RES_INTER = - 8 * MNI * MNJ * (MNI*MNI + MNJ*MNJ - 2 * MS*MS) / c / d * log((d + c) / (d - c));
    REAL RES = prefix * SQLamFactor * (RES_T + RES_U + RES_INTER);
    return RES;
}