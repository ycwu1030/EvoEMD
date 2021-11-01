#ifndef _BOLTZMANN_EQUATION_H_
#define _BOLTZMANN_EQUATION_H_

#include <string>
#include <vector>

#include "EvoEMD/Hubble.h"
#include "EvoEMD/ParticleBase.h"
#include "EvoEMD/RungeKutta.h"

namespace EvoEMD {

class BoltzmannEquation : public ODE_FUNCS {
private:
    Parameter_Base *ptr_scale;
    REAL scale;
    std::vector<int> poi_pids;
    std::vector<std::string> poi_names;
    std::vector<Pseudo_Particle *> poi_ptrs;
    Particle_Factory &pf;
    Hubble_History &hh;
    void Setup_Scale();

public:
    BoltzmannEquation(Parameter_Base *scale = nullptr);
    ~BoltzmannEquation(){};

    void Set_X_Range(REAL X_BEGIN, REAL X_END);
    virtual VD dYdX(REAL x, VD &y, VD &delta_y_ratio) override;
    virtual VD Yeq(REAL x) override;
    virtual VB Is_Thermalized() override;
    virtual VB Can_be_Negative() override;
};

}  // namespace EvoEMD

#endif  //_BOLTZMANN_EQUATION_H_
