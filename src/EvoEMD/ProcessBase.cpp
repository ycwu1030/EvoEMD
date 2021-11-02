#include "EvoEMD/ProcessBase.h"

#include "EvoEMD/spdlog_wrapper.h"

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

Process::~Process() {
    delete CR_Calculator;
    delete amp;
}

std::string Process::Get_Process_Name() const {
    std::string res = "";
    for (int i = 0; i < amp->N_INITIAL; i++) {
        res += amp->INITIAL[i]->Get_Name();
        res += " ";
    }
    res += "-> ";
    for (int i = 0; i < amp->N_FINAL; i++) {
        res += amp->FINAL[i]->Get_Name();
        res += " ";
    }
    return res;
}

REAL Process::Get_Collision_Rate(REAL T) { return CR_Calculator->Get_Collision_Rate(T); }

REAL Process::Get_Yield_Coeff(REAL T, int PID) { return amp->Get_Coeff(T, PID); }

Process_Factory::~Process_Factory() {
    for (auto &&proc : PL) {
        delete proc;
    }
}

Process_Factory &Process_Factory::Get_Process_Factory() {
    static Process_Factory pf;
    return pf;
}

}  // namespace EvoEMD