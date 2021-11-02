#include "ToyLG.h"

#include <cmath>
#include <iostream>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;
using namespace std;
int main(int argc, char const *argv[]) {
    // * Free Parameter can be modified accordingly
    // * Any other parameters that depend on these free parameters will update their value when their value is acquired
    Parameter_Base *ti = RETRIVE_PARAMETER(Ti);
    Parameter_Base *tr = RETRIVE_PARAMETER(Tr);
    Parameter_Base *mn1 = RETRIVE_PARAMETER(MN1);
    ti->Set_Value(15);
    tr->Set_Value(10);

    // *
    Pseudo_Particle *pp = RETRIVE_PARTICLE(900001);
    pp->Set_Init_Thermal_Status(false);
    BoltzmannEquation BE(mn1);
    REAL scale = mn1->Get_Value();
    REAL T_BEGIN = 100 * scale;
    REAL T_END = scale / 1000.0;
    BE.Set_X_Range(log(scale / T_BEGIN), log(scale / T_END));
    cout << "Solving for [" << BE.Get_X_BEGIN() << "," << BE.Get_X_END() << "]" << endl;
    cout << "DOF = " << BE.Get_DOF() << endl;
    cout << "Starting from: Y = " << BE.Get_Y_BEGIN() << endl;
    RungeKutta rkFI(&BE);
    // cout << "System Built" << endl;
    rkFI.Solve(1e-3, 1e-3);
    rkFI.Dump_Solution("ToyLG_FI_Result.txt");
    return 0;
}
