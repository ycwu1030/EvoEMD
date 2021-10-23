#include <iostream>

#include "Models/SeesawFreezeIn/Parameters.h"
using namespace EvoEMD;
using namespace std;
int main(int argc, char const *argv[]) {
    Parameter_Factory &pf = Parameter_Factory::Get_Parameter_Factory();
    double mass = pf.Get_Parameter_Value("MN1");
    cout << "MN1 mass = " << mass << endl;
    pf.Set_Parameter("MN1", 2000);
    mass = pf.Get_Parameter_Value("MN1");
    cout << "MN1 mass = " << mass << endl;
    return 0;
}
