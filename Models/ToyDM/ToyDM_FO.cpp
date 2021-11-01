#include "ToyDM_FO.h"

#include <cmath>
#include <iostream>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;
using namespace std;
int main(int argc, char const *argv[]) {
    Parameter_Base *tr = RETRIVE_PARAMETER(Tr);
    tr->Set_Value(1e6);
    Parameter_Base *mx = RETRIVE_PARAMETER(MX);
    BoltzmannEquation BE(mx);
    REAL scale = mx->Get_Value();
    REAL T_BEGIN = scale;
    REAL T_END = scale / 1000.0;
    BE.Set_X_Range(log(scale / T_BEGIN), log(scale / T_END));
    // BE.Set_X_BEGIN(log(scale / T_BEGIN) + 1e-6);
    cout << "Solving for [" << BE.Get_X_BEGIN() << "," << BE.Get_X_END() << "]" << endl;
    cout << "DOF = " << BE.Get_DOF() << endl;
    cout << "Starting from: Y = " << BE.Get_Y_BEGIN()[0] << endl;
    RungeKutta rk(&BE);
    cout << "System Built" << endl;
    rk.Solve(1e-3, 1e-3);
    rk.Dump_Solution("ToyDM_FO_Result.txt");
    return 0;
}
