#ifndef _PROCESS_BASE_H_
#define _PROCESS_BASE_H_

#include <set>

#include "EvoEMD/Cache.h"
#include "EvoEMD/CollisionRate.h"

namespace EvoEMD {
class Process {
public:
    using INITIAL_STATES = Amplitude::INITIAL_STATES;
    using FINAL_STATES = Amplitude::FINAL_STATES;

    static int NProcess;
    Process(Amplitude *amp);
    ~Process();

    std::string Get_Process_Name() const;
    int Get_Process_ID() const { return process_id; }
    REAL Get_Collision_Rate(REAL T);
    REAL Get_Yield_Coeff(REAL T, int PID);
    void Clean_Cache() { cr_cache.Clean_Cache(); }

protected:
    INITIAL_STATES INIT;
    FINAL_STATES FINAL;
    int process_id;
    Amplitude *amp;
    Collision_Rate *CR_Calculator;
    CACHE cr_cache;
};

class Process_Factory {
    // * This class is used to manage all process, it is a singleton
    // * The only object of this class will own all processes
    // * Users won't directly use this class unless one wants to dig into the process

public:
    typedef std::set<Process *> Process_List;
    static Process_Factory &Get_Process_Factory();
    static void Clean_Cache();
    void Register_Process(Process *proc) { PL.insert(proc); };
    Process_List Get_Process_List() { return PL; }

private:
    Process_List PL;

    Process_Factory(){};
    ~Process_Factory();
};

}  // namespace EvoEMD

class Register_Process {
public:
    Register_Process(EvoEMD::Process *proc) {
        EvoEMD::Process_Factory &pf = EvoEMD::Process_Factory::Get_Process_Factory();
        pf.Register_Process(proc);
    }
};

#define REGISTER_PROCESS(AmpClassName) \
    Register_Process g_register_process_##AmpClassName(new EvoEMD::Process(new AmpClassName))

#endif  //_PROCESS_BASE_H_
