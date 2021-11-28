#ifndef _TOY_DM_PROCESSES_H_
#define _TOY_DM_PROCESSES_H_

#include "EvoEMD/EvoEMD.h"

class XX_HH_Amp : public EvoEMD::Amplitude_Base {
private:
    REAL Sub1;

public:
    XX_HH_Amp();
    ~XX_HH_Amp(){};

    virtual void Update_Value(REAL input) override;
    virtual void Update_Amp(REAL sqrt_shat) override;
    virtual REAL Get_Offset(REAL T, int PID) override;
};

#endif  //_TOY_DM_PROCESSES_H_
