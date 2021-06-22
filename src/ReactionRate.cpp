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
    REAL gamma = LamX * LamX / 16.0 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) /
                 z * (1 + B / A - C / A) * SQLam;
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
    REAL A = pow(MN / MNR1, 2);
    REAL coup = real(Nu_Param.Get_YdagYij(i, i));
    REAL gamma = coup / 8 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) / z;
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
    // REAL s = smin + (smax - smin) * x[0];
    // REAL JacS = (smax - smin);
    REAL s = smin * exp(x[0] * log(smax / smin));
    REAL JacS = smin * log(smax / smin) * exp(x[0] * log(smax / smin));

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
    REAL prefix = 8 * M_PI * LamX * LamX * (s + MCHI * MCHI - MS * MS);
    REAL YY11 = real(Nu_Param.Get_YdagYij(1, 1));
    REAL YY22 = real(Nu_Param.Get_YdagYij(2, 2));
    REAL YY12 = real(Nu_Param.Get_YdagYij(1, 2));
    REAL N1Diag = (s + MNR1 * MNR1) / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2));
    REAL N2Diag = (s + MNR2 * MNR2) / (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterPref = 2 * (s + MNR1 * MNR2) / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2)) /
                     (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterR = InterPref * ((s - MNR1 * MNR1) * (s - MNR2 * MNR2) + MNR1 * MNR2 * GammaN1 * GammaN2);
    REAL RES = SQLamFactor * prefix * (YY11 * N1Diag + YY22 * N2Diag + YY12 * InterR);
    return RES;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_Kallen_NNChiChi(REAL s, int i, int j) {
    REAL MNI = i == 0 ? MNR1 : MNR2;
    REAL MNJ = j == 0 ? MNR1 : MNR2;
    if (s < pow(MNI + MNJ, 2) || s < pow(2 * MCHI, 2)) {
        return 0;
    }
    REAL SYM = i == j ? 1.0 / 2.0 : 1.0;
    REAL SQLamFactor =
        sqrt(Kallen_Lam(1.0, MNI * MNI / s, MNJ * MNJ / s) * Kallen_Lam(1.0, MCHI * MCHI / s, MCHI * MCHI / s));
    REAL prefix = 2 * M_PI * pow(LamX, 4) * SYM;
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
    REAL SYM = i == j ? 1.0 / 2.0 : 1.0;
    REAL SQLamFactor = sqrt(Kallen_Lam(1.0, MNI * MNI / s, MNJ * MNJ / s) * Kallen_Lam(1.0, MS * MS / s, MS * MS / s));
    REAL prefix = 2 * M_PI * pow(LamX, 4) * SYM;
    REAL EI = Ei(sqrt(s), MNI, MNJ);
    REAL EJ = Ei(sqrt(s), MNJ, MNI);
    REAL Pc = Pi(sqrt(s), MNI, MNJ);
    REAL E1 = Ei(sqrt(s), MS, MS);
    REAL E2 = Ei(sqrt(s), MS, MS);
    REAL Pf = Pi(sqrt(s), MS, MS);

    REAL c = 2 * Pc * Pf;
    REAL at = 2 * E1 * EI - MNI * MNI - MS * MS;
    REAL au = 2 * E1 * EJ - MNJ * MNJ - MS * MS;
    REAL bt = at + MCHI * MCHI;
    REAL bu = au + MCHI * MCHI;
    REAL k1 = MNI * MNI + MNJ * MNJ;
    REAL k22 = (MS * MS - MNI * MNI) * (MS * MS - MNJ * MNJ);
    // * T-Channel
    // * Int[((at-c*x)(au+c*x)-k1*(at-c*x)-k22)/(bt-c*x)^2,{x,-1,1}]
    REAL RES_T = -2 + (2 * bt + au - at - k1) / c * log((bt + c) / (bt - c)) +
                 2 * ((at - bt) * (au + bt - k1) - k22) / (bt * bt - c * c);
    // * U-Channel
    // * Int[((at-c*x)(au+c*x)-k1*(au+c*x)-k22)/(bu+c*x)^2,{x,-1,1}]
    REAL RES_U = -2 + (2 * bu + at - au - k1) / c * log((bu + c) / (bu - c)) +
                 2 * ((au - bu) * (at + bu - k1) - k22) / (bu * bu - c * c);
    // * Interference between t and u channel
    // * Int[1/(bt-c*x)/(bu+c*x),{x,-1,1}]
    REAL RES_INTER = 2 * MNI * MNJ * (MS * MS - MNI * MNI - MNJ * MNJ) / c / (bt + bu) *
                     (log((bt + c) / (bt - c)) + log((bu + c) / (bu - c)));
    return prefix * SQLamFactor * (RES_T + RES_U + RES_INTER);
}
