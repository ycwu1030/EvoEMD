#include "ToyLG.h"

#include <cmath>
#include <iostream>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;
using namespace std;
int main(int argc, char const *argv[]) {
    Parameter_Base *ti = RETRIVE_PARAMETER(Ti);
    ti->Set_Value(15);
    Parameter_Base *mn1 = RETRIVE_PARAMETER(MN1);
    BoltzmannEquation BE(mn1);
    REAL scale = mn1->Get_Value();
    REAL T_BEGIN = 100 * scale;
    REAL T_END = scale / 1000.0;
    BE.Set_X_Range(log(scale / T_BEGIN), log(scale / T_END));
    cout << "Solving for [" << BE.Get_X_BEGIN() << "," << BE.Get_X_END() << "]" << endl;
    cout << "DOF = " << BE.Get_DOF() << endl;
    cout << "Starting from: Y = " << BE.Get_Y_BEGIN() << endl;
    RungeKutta rk(&BE);
    // cout << "System Built" << endl;
    rk.Solve(1e-3, 1e-3);
    rk.Dump_Solution("ToyLG_Result.txt");
    return 0;
}
