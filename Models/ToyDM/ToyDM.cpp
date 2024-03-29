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
    Parameter_Base *hid = RETRIEVE_PARAMETER(HubbleMethod);

    // * Ti and Tr can be reset to any value
    ti->Set_Value(1e5);
    tr->Set_Value(1);
    // * Any particle can also be accessed
    Particle_Base *pp = RETRIEVE_PARTICLE(900001);

    Boltzmann_Equation BE(mx);
    REAL scale = mx->Get_Value();

    // * Using default `Splitting' method for Hubble parameter
    // * For Freeze-Out
    REAL T_BEGIN = scale;
    REAL T_END = 1e-2;
    BE.Set_T_Range(T_BEGIN, T_END);
    pp->Set_Init_Thermal_Status(true);
    BE.Solve(1e-3, 1e-3);
    BE.Dump_Solution("ToyDM_FO_Result_SP.txt");

    // * For Freeze-In
    lam->Set_Value(1e-10);
    scale = mx->Get_Value();
    T_BEGIN = 100 * scale;
    T_END = 1e-2;
    BE.Set_T_Range(T_BEGIN, T_END);
    pp->Set_Init_Thermal_Status(false);
    BE.Solve(1e-3, 1e-5);
    BE.Dump_Solution("ToyDM_FI_Result_SP.txt");

    // * Using `BE' method for Hubble parameter
    hid->Set_Value(1);

    // * For Freeze-Out
    lam->Set_Value(0.4);
    T_BEGIN = scale;
    T_END = 1e-2;
    BE.Set_T_Range(T_BEGIN, T_END);
    pp->Set_Init_Thermal_Status(true);
    BE.Solve(1e-3, 1e-3);
    BE.Dump_Solution("ToyDM_FO_Result_BE.txt");

    // * For Freeze-In
    lam->Set_Value(1e-10);
    scale = mx->Get_Value();
    T_BEGIN = 100 * scale;
    T_END = 1e-2;
    BE.Set_T_Range(T_BEGIN, T_END);
    pp->Set_Init_Thermal_Status(false);
    BE.Solve(1e-3, 1e-5);
    BE.Dump_Solution("ToyDM_FI_Result_BE.txt");
    return 0;
}
