#include "EvoEMD/BoltzmannEquation.h"

#include <cmath>

namespace EvoEMD {

BoltzmannEquation::BoltzmannEquation(Parameter_Base *ptr_scale_)
    : pf(EvoEMD::Particle_Factory::Get_Particle_Factory()), ptr_scale(ptr_scale_) {
    std::set<int> POI = pf.Get_POI();
    DOF = POI.size();
    std::set<int>::iterator iter = POI.begin();
    for (; iter != POI.end(); iter++) {
        poi_pids.push_back(*iter);
        Pseudo_Particle *pp = pf.Get_Particle(*iter);
        poi_ptrs.push_back(pp);
        poi_names.push_back(pp->Get_Name());
    }
    Setup_Scale();
}

void BoltzmannEquation::Setup_Scale() {
    if (!ptr_scale) {
        scale = 1000.0;
    } else {
        scale = ptr_scale->Get_Value();
    }
}

VD BoltzmannEquation::dYdX(REAL x, VD y) {
    // * x = log(z), z = scale/T
    // * Then the Boltzmann equation for particle Yield evolution is
    // * T^3H(beta_T dY/dx + 3(1-beta_T)Y) = CollisionRate
    // * beta_T indicates the dependence of T on scale factor: T~ a^{-beta_T}
    REAL z = exp(x);
    REAL T = scale / z;
}

}  // namespace EvoEMD
