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
    REAL T_BEGIN = scale / 40.0;
    REAL T_END = scale / 100.0;
    BE.Set_X_BEGIN(log(scale / T_BEGIN));
    BE.Set_X_END(log(scale / T_END));
    VD BD(1);
    BD[0] = (Particle_Factory::Get_Particle_Factory().Get_Particle(900001)->Get_Equilibrium_Yield_at_T(T_BEGIN));
    BE.Set_BOUNDARY_CONDITION(BD);
    // BE.Set_X_BEGIN(log(scale / T_BEGIN) + 1e-6);
    BE.dYdX(log(scale / T_BEGIN), BD);
    cout << "Solving for [" << BE.Get_X_BEGIN() << "," << BE.Get_X_END() << "]" << endl;
    cout << "DOF = " << BE.Get_DOF() << endl;
    cout << "Starting from: Y = " << BE.Get_BOUNDARY_CONDITION()[0] << endl;
    RungeKutta rk(&BE);
    cout << "System Built" << endl;
    rk.Solve(1e-2, 1e-3);
    rk.Dump_Solution("ToyDM_Result.txt");
    return 0;
}
