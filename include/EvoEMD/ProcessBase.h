#ifndef _PROCESS_BASE_H_
#define _PROCESS_BASE_H_

#include "EvoEMD/CollisionRate.h"

namespace EvoEMD {
class Process {
public:
    using INITIAL_STATES = Amplitude::INITIAL_STATES;
    using FINAL_STATES = Amplitude::FINAL_STATES;

    Process(Amplitude *amp);
    ~Process();

    REAL Get_Collision_Rate(REAL T);
    REAL Get_Yield_Coeff(REAL T);

protected:
    INITIAL_STATES INIT;
    FINAL_STATES FINAL;
    Amplitude *amp;
    Collision_Rate *CR_Calculator;
};
}  // namespace EvoEMD
#endif  //_PROCESS_BASE_H_
