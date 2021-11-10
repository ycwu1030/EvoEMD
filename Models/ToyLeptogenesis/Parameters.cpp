#include "Parameters.h"

param_GammaN1::param_GammaN1() : Parameter_Base("GammaN1") {
    Parameter_Base* p_mn1 = RETRIVE_PARAMETER(MN1);
    Parameter_Base* p_lam = RETRIVE_PARAMETER(Lam);
    parent_parameters.insert(p_mn1);
    parent_parameters.insert(p_lam);
    p_mn1->Register_Descendent_Parameter(this);
    p_lam->Register_Descendent_Parameter(this);
}

void param_GammaN1::Update_Value(REAL input) {
    REAL mn1 = GET_PARAMETER_VALUE(MN1);
    REAL lam = GET_PARAMETER_VALUE(Lam);
    value = lam * lam * mn1 / 8 / M_PI;
}
