#include "EvoEMD/BoltzmannEquation.h"

namespace EvoEMD {

BoltzmannEquation::BoltzmannEquation() : pf(EvoEMD::Particle_Factory::Get_Particle_Factory()) {
    std::set<int> POI = pf.Get_POI();
    DOF = POI.size();
    std::set<int>::iterator iter = POI.begin();
    for (; iter != POI.end(); iter++) {
        poi_pids.push_back(*iter);
        Pseudo_Particle *pp = pf.Get_Particle(*iter);
        poi_ptrs.push_back(pp);
    }
}

}  // namespace EvoEMD
