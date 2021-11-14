#ifndef _TOY_LEPTOGENESIS_AMP_H_
#define _TOY_LEPTOGENESIS_AMP_H_

#include "EvoEMD/EvoEMD.h"

class N_LPhi_Amp : public EvoEMD::Amplitude_Base {
public:
    N_LPhi_Amp();

    virtual void Update_Amp(REAL sqrt_shat) override;
    virtual REAL Get_Coeff(REAL T, int PID) override;
};

class delta_N_LPhi_Amp : public EvoEMD::Amplitude_Base {
public:
    delta_N_LPhi_Amp();
    virtual void Update_Amp(REAL sqrt_shat) override;
    virtual REAL Get_Coeff(REAL T, int PID) override;
};

#endif  //_TOY_LEPTOGENESIS_AMP_H_
