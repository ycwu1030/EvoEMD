#ifndef _TOY_LEPTOGENESIS_AMP_H_
#define _TOY_LEPTOGENESIS_AMP_H_

#include "EvoEMD/EvoEMD.h"

class N_LPhi_Amp_CPC : public EvoEMD::Amplitude_Base {
private:
    REAL Sub1;

public:
    N_LPhi_Amp_CPC();

    virtual void Update_Value(REAL input) override;
    virtual void Update_Amp(REAL sqrt_shat) override;
    virtual REAL Get_Offset(REAL T, int PID) override;
};

class N_LPhi_Amp_CPV : public EvoEMD::Amplitude_Base {
private:
    REAL Sub1;

public:
    N_LPhi_Amp_CPV();

    virtual void Update_Value(REAL input) override;
    virtual void Update_Amp(REAL sqrt_shat) override;
    virtual REAL Get_Offset(REAL T, int PID) override;
};

#endif  //_TOY_LEPTOGENESIS_AMP_H_
