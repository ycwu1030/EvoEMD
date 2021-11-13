#include "ToyLG.h"

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
    Parameter_Base *mn1 = RETRIEVE_PARAMETER(MN1);
    ti->Set_Value(15);
    tr->Set_Value(10);

    // * Any particle can be accessed
    Particle_Base *pp = RETRIEVE_PARTICLE(900001);

    auto pl = Process_Factory::Get_Process_Factory().Get_Process_List();
    for (auto &&pc : pl) {
        cout << "PID: " << pc->Get_Process_ID() << " " << pc->Get_Process_Name() << endl;
    }

    Boltzmann_Equation BE(mn1);
    REAL scale = mn1->Get_Value();
    REAL T_BEGIN = 100 * scale;
    REAL T_END = scale / 1000.0;
    BE.Set_T_Range(T_BEGIN, T_END);
    pp->Set_Init_Thermal_Status(true);
    BE.Solve(1e-3, 1e-3);
    BE.Dump_Solution("ToyLG_FO_Result.txt");

    pp->Set_Init_Thermal_Status(false);
    BE.Solve(1e-3, 1e-3);
    BE.Dump_Solution("ToyLG_FI_Result.txt");

    return 0;
}
