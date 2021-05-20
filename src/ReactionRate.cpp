#include "LeptogenesisRate.h"
#include "cuba.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_sf_bessel.h"
#include "spdlog_wrapper.h"
#include <cmath>
#include <algorithm>

using namespace std;

// CUBA PARAMETERS:
// COMMON VARIABLES
#define MINEVAL 10000
#define MAXEVAL 100000
#define EPSREL 1E-3
#define EPSABS 0 // 1E-7 precision is about 1 ab, which is way small.
#define NVEC 1
#define FLAGS 0
#define SEED 0

// VEGAS
#define NSTART 2000
#define NINCREASE 500
#define NBATCH 2000
#define GRIDNO 0

REAL LeptogenesisRate::Calc_NChiS_Gamma(REAL Temp, int i)
{
    REAL MN = i == 0 ? MNR1 : MNR2;
    if (MN < MCHI + MS)
    {
        return 0;
    }
    REAL A = MN / MNR1;
    REAL B = MCHI / MNR1;
    REAL C = MS / MNR1;
    REAL z = MNR1 / Temp;
    REAL SQLam = SQkaller(1, B / A, C / A);
    REAL gamma = LamX * LamX / 16 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) / z * (1 + B / A - C / A) * SQLam;
    return gamma;
}

REAL LeptogenesisRate::Calc_ChiSN_Gamma(REAL Temp, int i)
{
    REAL MN = i == 0 ? MNR1 : MNR2;
    if (MCHI < MN + MS)
    {
        return 0;
    }
    REAL A = MCHI / MNR1;
    REAL B = MN / MNR1;
    REAL C = MS / MNR1;
    REAL z = MNR1 / Temp;
    REAL SQLam = SQkaller(1, B / A, C / A);
    REAL gamma = LamX * LamX / 16 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) / z * (1 + B / A - C / A) * SQLam;
    return gamma;
}

REAL LeptogenesisRate::Calc_SChiN_Gamma(REAL Temp, int i)
{
    REAL MN = i == 0 ? MNR1 : MNR2;
    if (MS < MN + MCHI)
    {
        return 0;
    }
    REAL A = MS / MNR1;
    REAL B = MN / MNR1;
    REAL C = MCHI / MNR1;
    REAL z = MNR1 / Temp;
    REAL SQLam = SQkaller(1, B / A, C / A);
    REAL gamma = LamX * LamX / 16 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) / z * (1 - B / A - C / A) * SQLam;
    return gamma;
}

REAL LeptogenesisRate::Calc_NLPhi_Gamma(REAL Temp, int i)
{
    REAL z = MNR1 / Temp;
    REAL MN = i == 0 ? MNR1 : MNR2;
    REAL A = MN / MNR1;
    REAL coup = real(Nu_Param.Get_YdagYij(i, i));
    REAL gamma = coup * coup / 8 / pow(M_PI, 3) * pow(MNR1, 4) * sqrt(pow(A, 3)) * gsl_sf_bessel_K1(z * sqrt(A)) / z;
    return 1;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_SQkaller(REAL s, int processid, int i)
{
    switch (processid)
    {
    case 0:
        return SqAmp_dOmega_with_SQkaller_LPhiChiS(s);
    case 1:
        return SqAmp_dOmega_with_SQkaller_NNChiChi(s, i);
    case 2:
        return SqAmp_dOmega_with_SQkaller_NNSS(s, i);
    default:
        return SqAmp_dOmega_with_SQkaller_LPhiChiS(s);
    }
}

struct INTE_PARAM
{
    LeptogenesisRate *ptr;
    int processid;
    int nrid;
    REAL Temp;
    REAL smin;
    REAL smax;
};

int gamma_Integrand(const int *ndim, const REAL x[], const int *ncomp, REAL ff[], void *_params)
{
    INTE_PARAM *params = (INTE_PARAM *)_params;
    int procid = params->processid;
    int nrid = params->nrid;
    REAL Temp = params->Temp;
    REAL smin = params->smin;
    REAL smax = params->smax;
    REAL s = smin + (smax - smin) * x[0];
    ff[0] = (smax - smin) * params->ptr->SqAmp_dOmega_with_SQkaller(s, procid, nrid) * sqrt(s) * gsl_sf_bessel_K1(sqrt(s) / Temp) * Temp / 32.0 / pow(2 * M_PI, 6);
    return 0;
}

REAL LeptogenesisRate::Calc_LPhiChiS_Gamma(REAL Temp)
{
    INTE_PARAM param = {this, 0, -1, Temp, pow(MCHI + MS, 2), pow(100 * Temp, 2)};
    int fail, neval;
    const int NDIM = 1;
    const int NCOMP = 1;
    REAL RES[NCOMP];
    REAL ERR[NCOMP];
    REAL PROB[NCOMP];

    Vegas(NDIM, NCOMP, gamma_Integrand, &param, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);

    return RES[0];
}

REAL LeptogenesisRate::Calc_NNChiChi_Gamma(REAL Temp, int i)
{
    REAL MN = i == 0 ? MNR1 : MNR2;
    REAL mmin = max(2 * MN, 2 * MCHI);
    INTE_PARAM param = {this, 1, i, Temp, pow(mmin, 2), pow(100 * Temp, 2)};
    int fail, neval;
    const int NDIM = 1;
    const int NCOMP = 1;
    REAL RES[NCOMP];
    REAL ERR[NCOMP];
    REAL PROB[NCOMP];

    Vegas(NDIM, NCOMP, gamma_Integrand, &param, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);

    return RES[0];
}

REAL LeptogenesisRate::Calc_NNSS_Gamma(REAL Temp, int i)
{
    REAL MN = i == 0 ? MNR1 : MNR2;
    REAL mmin = max(2 * MN, 2 * MCHI);
    INTE_PARAM param = {this, 2, i, Temp, pow(mmin, 2), pow(100 * Temp, 2)};
    int fail, neval;
    const int NDIM = 1;
    const int NCOMP = 1;
    REAL RES[NCOMP];
    REAL ERR[NCOMP];
    REAL PROB[NCOMP];

    Vegas(NDIM, NCOMP, gamma_Integrand, &param, NVEC, EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, GRIDNO, NULL, NULL, &neval, &fail, RES, ERR, PROB);

    return RES[0];
}

REAL LeptogenesisRate::SqAmp_dOmega_with_SQkaller_LPhiChiS(REAL s)
{
    REAL M1 = 0;
    REAL M2 = 0;
    REAL M3 = MCHI;
    REAL M4 = MS;
    if (s < pow(M3 + M4, 2) || s < pow(M1 + M2, 2))
    {
        return 0;
    }
    REAL SQLamFactor = sqrt(SQkaller(1.0, M1 * M1 / s, M2 * M2 / s) * SQkaller(1.0, M3 * M3 / s, M3 * M3 / s));
    REAL prefix = 2 * M_PI * LamX * LamX * (s + MCHI * MCHI - MS * MS);
    REAL h1a2coup = 0;
    REAL h2a2coup = 0;
    REAL Rh1aha2coup = 0;
    REAL Ih1aha2coup = 0;
    for (int lep_index = 0; lep_index < 3; lep_index++)
    {
        auto h1a = Nu_Param.Get_Yij(0, lep_index);
        auto h2a = Nu_Param.Get_Yij(1, lep_index);
        h1a2coup += norm(h1a); // norm is the square of its magnitude
        h2a2coup += norm(h2a);
        Rh1aha2coup += real(h1a * conj(h2a));
        Ih1aha2coup += imag(h1a * conj(h2a));
    }
    REAL N1Diag = (s + MNR1 * MNR1) / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2));
    REAL N2Diag = (s + MNR2 * MNR2) / (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterPref = 2 * (s + MNR1 * MNR2) / (pow(s - MNR1 * MNR1, 2) + pow(MNR1 * GammaN1, 2)) / (pow(s - MNR2 * MNR2, 2) + pow(MNR2 * GammaN2, 2));
    REAL InterR = InterPref * ((s - MNR1 * MNR1) * (s - MNR2 * MNR2) + MNR1 * MNR2 * GammaN1 * GammaN2);
    REAL InterI = InterPref * ((s - MNR1 * MNR1) * MNR2 * GammaN2 - (s - MNR2 * MNR2) * MNR1 * GammaN1);
    REAL RES = SQLamFactor * prefix * (h1a2coup * N1Diag + h2a2coup * N2Diag + Rh1aha2coup * InterR + Ih1aha2coup * InterI);
    return RES;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_SQkaller_NNChiChi(REAL s, int i)
{
    return 0;
}

REAL LeptogenesisRate::SqAmp_dOmega_with_SQkaller_NNSS(REAL s, int i)
{
    return 0;
}
