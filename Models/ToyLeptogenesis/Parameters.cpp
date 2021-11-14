#include "Parameters.h"

param_GammaN1::param_GammaN1() : Parameter_Base("GammaN1") {
    Parameter_Base* p_mn1 = RETRIEVE_PARAMETER(MN1);
    Parameter_Base* p_lam = RETRIEVE_PARAMETER(Lam);

    Register_Dependencies(p_mn1, p_lam);
}

void param_GammaN1::Update_Value(REAL input) {
    REAL mn1 = GET_PARAMETER_VALUE(MN1);
    REAL lam = GET_PARAMETER_VALUE(Lam);
    value = lam * lam * mn1 / 8 / M_PI;
}
