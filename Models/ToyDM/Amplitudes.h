#ifndef _TOY_DM_PROCESSES_H_
#define _TOY_DM_PROCESSES_H_

#include "EvoEMD/EvoEMD.h"

class XX_SS_Amp : public EvoEMD::Amplitude_Base {
private:
    REAL Sub1;

public:
    XX_SS_Amp();
    ~XX_SS_Amp(){};

    virtual void Update_Value(REAL input) override;
    virtual void Update_Amp(REAL sqrt_shat) override;
    virtual REAL Get_Coeff(REAL T, int PID) override;
};

#endif  //_TOY_DM_PROCESSES_H_
