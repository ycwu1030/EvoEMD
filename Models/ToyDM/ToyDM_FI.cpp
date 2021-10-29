#include "ToyDM_FI.h"

#include <cmath>
#include <iostream>

#include "EvoEMD/EvoEMD.h"

using namespace EvoEMD;
using namespace std;
int main(int argc, char const *argv[]) {
    Parameter_Base *tr = RETRIVE_PARAMETER(Tr);
    tr->Set_Value(1);
    REAL Ti = GET_PARAM_VALUE(Ti);
    Parameter_Base *mx = RETRIVE_PARAMETER(MX);
    BoltzmannEquation BE(mx);
    REAL scale = mx->Get_Value();
    REAL T_BEGIN = 100 * scale;
    REAL T_END = 1e-2;
    BE.Set_X_Range(log(scale / T_BEGIN), log(scale / T_END));
    // BE.Set_X_BEGIN(log(scale / T_BEGIN) + 1e-6);
    cout << "Solving for [" << BE.Get_X_BEGIN() << "," << BE.Get_X_END() << "]" << endl;
    cout << "DOF = " << BE.Get_DOF() << endl;
    cout << "Starting from: Y = " << BE.Get_Y_BEGIN()[0] << endl;
    Hubble_History::Get_Hubble_History().Print_History();
    RungeKutta rk(&BE);
    rk.Solve(1e-4);
    rk.Dump_Solution("ToyDM_FI_Result.txt");
    return 0;
}
