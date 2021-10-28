#include <cmath>
#include <iostream>

#include "EvoEMD/BoltzmannEquation.h"
#include "EvoEMD/RungeKutta.h"
#include "Models/ToyDM/Parameters.h"
#include "Models/ToyDM/Particles.h"
#include "Models/ToyDM/Processes.h"
#include "gsl/gsl_sf_bessel.h"

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
    rk.Solve(1e-4);
    rk.Dump_Solution("ToyDM_Result.txt");
    return 0;
}
