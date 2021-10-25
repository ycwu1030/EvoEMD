#include "EvoEMD/BoltzmannEquation.h"

#include <cmath>

#include "EvoEMD/ProcessBase.h"

namespace EvoEMD {

BoltzmannEquation::BoltzmannEquation(Parameter_Base *ptr_scale_)
    : pf(EvoEMD::Particle_Factory::Get_Particle_Factory()),
      hh(EvoEMD::Hubble_History::Get_Hubble_History()),
      ptr_scale(ptr_scale_) {
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

    for (int i = 0; i < DOF; i++) {
        poi_ptrs[i]->Yield = y[i];
    }

    REAL z = exp(x);
    REAL T = scale / z;
    int hubble_period_id = hh.Get_Period_ID_at_T(T);
    Hubble_For_Single_Period *hs = hh[hubble_period_id];
    double beta_T = hs->Get_beta_T();
    REAL HatT = hs->Get_Hubble_at_T(T);
    bool isentropic = hs->Is_Isentropic();
    VD res(DOF, 0);
    for (int i = 0; i < DOF; i++) {
        Pseudo_Particle *pp = poi_ptrs[i];
        std::set<Process *> sp = pp->Get_Process();
        for (auto &&proc_ptr : sp) {
            res[i] += proc_ptr->Get_Collision_Rate(T) * proc_ptr->Get_Yield_Coeff(T, pp->Get_PID());
        }
    }
    for (int i = 0; i < DOF; i++) {
        res[i] *= 1.0 / pow(T, 3) / HatT;
        if (!isentropic) {
            // that beta_T != 1
            res[i] -= 3.0 * (1.0 - beta_T) * y[i];
            res[i] /= beta_T;
        }  // else beta_T == 1, no further action needed
    }
    return res;
}

}  // namespace EvoEMD
