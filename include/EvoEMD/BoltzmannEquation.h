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
    REAL T_BEGIN;
    REAL T_END;
    REAL scale;
    std::vector<int> poi_pids;
    std::vector<std::string> poi_names;
    std::vector<Pseudo_Particle *> poi_ptrs;
    Particle_Factory &pf;
    Hubble_History &hh;
    RungeKutta rk;
    void Setup_Scale();
    void Set_X_Range(REAL X_BEGIN, REAL X_END);
    REAL dYidX(int i, REAL x, const VD &y, const VD &delta_y_ratio);

public:
    BoltzmannEquation(Parameter_Base *scale = nullptr);
    ~BoltzmannEquation(){};

    // * User may not want to use following methods
    virtual VD dYdX(REAL x, const VD &y, const VD &delta_y_ratio) override;
    virtual VD Yeq(REAL x) override;
    virtual VB Is_Thermalized() override;
    virtual VB Can_be_Negative() override;
    virtual VB Should_be_Thermalized(REAL x, const VD &y, const VD &delta_y_ratio) override;

    // * User should only use following methods
    void Set_T_Range(REAL T_BEGIN, REAL T_END);
    void Set_Z_Range(REAL Z_BEGIN, REAL Z_END);
    RungeKutta::STATUS Solve(REAL step_size, REAL eps_rel = 1e-4);
    void Dump_Solution(std::string filename);
};

}  // namespace EvoEMD

#endif  //_BOLTZMANN_EQUATION_H_