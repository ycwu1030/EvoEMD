#ifndef _TOY_DM_PROCESSES_H
#define _TOY_DM_PROCESSES_H

#include "EvoEMD/EvoEMD.h"

class XX_SS_Amp : public EvoEMD::Amplitude {
public:
    XX_SS_Amp();
    ~XX_SS_Amp(){};

    virtual void Update_Amp(REAL sqrt_shat) override;
    virtual REAL Get_Coeff(REAL T, int PID) override;
};

#endif  //_TOY_DM_PROCESSES_H
