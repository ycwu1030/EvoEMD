#include "EvoEMD/ProcessBase.h"

namespace EvoEMD {

Process::Process(Amplitude *amp) {
    this->amp = amp;
    INIT = amp->INITIAL;
    FINAL = amp->FINAL;
    if (amp->N_INITIAL == 1 && amp->N_FINAL == 2) {
        CR_Calculator = new Decay12_Rate(amp);
    }
    if (amp->N_INITIAL == 2 && amp->N_FINAL == 2) {
        CR_Calculator = new Scatter22_Rate(amp);
    }
    for (int i = 0; i < INIT.size(); i++) {
        INIT[i]->Register_Process(this);
    }
    for (int i = 0; i < FINAL.size(); i++) {
        FINAL[i]->Register_Process(this);
    }
}

Process::~Process() { delete CR_Calculator; }

REAL Process::Get_Collision_Rate(REAL T) { return CR_Calculator->Get_Collision_Rate(T); }

REAL Process::Get_Yield_Coeff(REAL T, int PID) { return amp->Get_Coeff(T,PID); }

}  // namespace EvoEMD
