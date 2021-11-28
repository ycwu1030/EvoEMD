#include "EvoEMD/ProcessBase.h"

#include "EvoEMD/Cache.h"
#include "EvoEMD/spdlog_wrapper.h"

namespace EvoEMD {

int Process::NProcess = 0;

Process::Process(Amplitude_Base *amp) {
    this->amp = amp;
    INIT = amp->INITIAL;
    FINAL = amp->FINAL;
    if (amp->N_INITIAL == 1 && amp->N_FINAL == 2) {
        CR_Calculator = new Decay_Rate(amp);
    }
    if (amp->N_INITIAL == 2 && amp->N_FINAL == 2) {
        CR_Calculator = new Scatter_Rate(amp);
    }
    for (int i = 0; i < INIT.size(); i++) {
        INIT[i]->Register_Process(this);
    }
    for (int i = 0; i < FINAL.size(); i++) {
        FINAL[i]->Register_Process(this);
    }
    process_id = NProcess++;
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

REAL Process::Get_Collision_Rate(REAL T) {
    REAL res;
    INDEX ind = OBTAIN_KEY(T);
    bool exist = cr_cache.Get(ind, res);
    if (!exist) {
        res = CR_Calculator->Get_Collision_Rate(T);
        cr_cache.Insert(ind, res);
    }
    return res;
}

REAL Process::Get_Offset(REAL T, int PID) { return amp->Get_Offset(T, PID); }

Process_Factory::~Process_Factory() {
    for (auto &&proc : PL) {
        delete proc;
    }
}

Process_Factory &Process_Factory::Get_Process_Factory() {
    static Process_Factory pf;
    return pf;
}

void Process_Factory::Clean_Cache() {
    Process_List pl = Get_Process_Factory().Get_Process_List();
    for (auto &&prc : pl) {
        prc->Clean_Cache();
    }
}

}  // namespace EvoEMD
