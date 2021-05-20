#include "LeptogenesisRate.h"

REAL LeptogenesisRate::Calc_NChiS_Gamma(REAL Temp, REAL MA, REAL MB, REAL MC)
{
    if (MA <= MB + MC)
    {
        return 0;
    }
    return 1;
}

REAL LeptogenesisRate::Calc_NLPhi_Gamma(REAL Temp, int i)
{
    double MN = i == 0 ? MNR1 : MNR2;
    return 1;
}

REAL LeptogenesisRate::Calc_LPhiChiS_Gamma(REAL Temp)
{
    return 1;
}

REAL LeptogenesisRate::Calc_NNChiChi_Gamma(REAL Temp, int i)
{
    double MN = i == 0 ? MNR1 : MNR2;
    return 1;
}

REAL LeptogenesisRate::Calc_NNSS_Gamma(REAL Temp, int i)
{
    double MN = i == 0 ? MNR1 : MNR2;
    return 1;
}
