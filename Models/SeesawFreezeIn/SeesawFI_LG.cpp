#include <iostream>

#include "EvoEMD/Hubble.h"
#include "Models/SeesawFreezeIn/Parameters.h"
#include "Models/SeesawFreezeIn/Particles.h"
#include "Models/SeesawFreezeIn/SeesawTypeI.h"
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
    // pf.List_Parameters();
    Parameter_Base *pb = RETRIVE_PARAMETER(Seesaw_Parameters);
    cout << " YdagY(0,0) = " << ((SeesawTypeI *)pb)->Get_YdagYij(0, 0) << endl;
    pf.Set_Parameter_Value("MN1", 3001.1);
    cout << " YdagY(0,0) = " << ((SeesawTypeI *)pb)->Get_YdagYij(0, 0) << endl;
    Hubble_History &hh = Hubble_History::Get_Hubble_History();
    hh.Print_History();
    cout << " H(T=30000) = " << hh.Get_Hubble_at_T(30000) << endl;
    return 0;
}