#include "ToyDM.h"

#include <cmath>
#include <iostream>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;
using namespace std;
int main(int argc, char const *argv[]) {
    // * Free Parameter can be modified accordingly
    // * Any other parameters that depend on these free parameters will update their value when their value is acquired
    Parameter_Base *ti = RETRIEVE_PARAMETER(Ti);
    Parameter_Base *tr = RETRIEVE_PARAMETER(Tr);
    Parameter_Base *mx = RETRIEVE_PARAMETER(MX);
    Parameter_Base *lam = RETRIEVE_PARAMETER(Lam);

    // * Ti and Tr can be reset to any value
    ti->Set_Value(10);
    tr->Set_Value(1);

    // * Any particle can also be accessed
    Pseudo_Particle *pp = RETRIVE_PARTICLE(900001);

    Boltzmann_Equation BE(mx);
    REAL scale = mx->Get_Value();

    // * For Freeze-Out
    REAL T_BEGIN = scale;
    REAL T_END = scale / 1000.0;
    BE.Set_T_Range(T_BEGIN, T_END);
    pp->Set_Init_Thermal_Status(true);
    BE.Solve(1e-3, 1e-3);
    BE.Dump_Solution("ToyDM_FO_Result.txt");

    // * For Freeze-In
    lam->Set_Value(1e-10);
    scale = mx->Get_Value();
    T_BEGIN = 100 * scale;
    T_END = 1e-2;
    BE.Set_T_Range(T_BEGIN, T_END);
    pp->Set_Init_Thermal_Status(false);
    BE.Solve(1e-3, 1e-3);
    BE.Dump_Solution("ToyDM_FI_Result.txt");
    return 0;
}
