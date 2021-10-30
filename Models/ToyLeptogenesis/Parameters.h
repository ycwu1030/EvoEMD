#ifndef _TOY_LEPTOGENESIS_PARAMETER_H_
#define _TOY_LEPTOGENESIS_PARAMETER_H_

#include "EvoEMD/EvoEMD.h"
using namespace EvoEMD;

class param_GammaN1 : public Parameter_Base {
public:
    param_GammaN1();

    virtual void Update_Value(REAL input) override;
};

#endif  //_TOY_LEPTOGENESIS_PARAMETER_H_
