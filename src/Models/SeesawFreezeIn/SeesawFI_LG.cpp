#include <iostream>

#include "Models/SeesawFreezeIn/Parameters.h"
#include "Models/SeesawFreezeIn/Particles.h"
using namespace EvoEMD;
using namespace std;
int main(int argc, char const *argv[]) {
    Parameter_Factory &pf = Parameter_Factory::Get_Parameter_Factory();
    double mass = pf.Get_Parameter_Value("MN1");
    cout << "MN1 mass = " << mass << endl;
    pf.Set_Parameter_Value("MN1", 2000);
    mass = pf.Get_Parameter_Value("MN1");
    cout << "MN1 mass = " << mass << endl;
    Particle_Factory &partF = Particle_Factory::Get_Particle_Factory();
    mass = partF.Get_Particle(900001)->Get_Mass();
    cout << "Particle N1 with PID = 900001 has mass " << mass << endl;
    return 0;
}
