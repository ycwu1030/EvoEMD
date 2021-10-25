#ifndef _BOLTZMANN_EQUATION_H_
#define _BOLTZMANN_EQUATION_H_

#include <string>
#include <vector>

#include "EvoEMD/ParticleBase.h"
#include "EvoEMD/RungeKutta.h"

namespace EvoEMD {

class BoltzmannEquation : public ODE_FUNCS {
private:
    std::vector<int> poi_pids;
    std::vector<std::string> poi_names;
    std::vector<Pseudo_Particle *> poi_ptrs;
    Particle_Factory &pf;

public:
    BoltzmannEquation();
    ~BoltzmannEquation();

    virtual VD dYdX(REAL x, VD y) override;
};

}  // namespace EvoEMD

#endif  //_BOLTZMANN_EQUATION_H_
