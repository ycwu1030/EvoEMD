#ifndef _BOLTZMANN_EQUATION_H_
#define _BOLTZMANN_EQUATION_H_

#include <string>
#include <vector>

#include "EvoEMD/Hubble.h"
#include "EvoEMD/ParticleBase.h"
#include "EvoEMD/RungeKutta.h"

namespace EvoEMD {

class Boltzmann_Equation : public ODE_FUNCS {
private:
    Parameter_Base *ptr_scale;
    REAL T_BEGIN;
    REAL T_END;
    REAL scale;
    std::vector<int> poi_pids;
    std::vector<std::string> poi_names;
    std::vector<Particle_Base *> poi_ptrs;
    Particle_Factory &pf;
    Hubble_History &hh;
    RungeKutta rk;
    void Setup_Scale();
    void Set_X_Range(REAL X_BEGIN, REAL X_END);
    REAL dYidX(int i, REAL x, const VD &y, const VD &delta_y_ratio);

public:
    /**
     * @brief default ctor for Boltzmann equation
     * @note   The value of the scale is used to compute z = m/T
     * @param  *scale: a ptr to the scale parameter
     * @retval
     */
    Boltzmann_Equation(Parameter_Base *scale = nullptr);
    ~Boltzmann_Equation(){};

    // * User may not want to use following methods
    virtual VD dYdX(REAL x, const VD &y, const VD &delta_y_ratio) override;
    virtual VD Yeq(REAL x) override;
    virtual VB Is_Thermalized() override;
    virtual VB Can_be_Negative() override;
    virtual VB Should_be_Thermalized(REAL x, const VD &y, const VD &delta_y_ratio) override;

    // * User should only use following methods
    /**
     * @brief  Setting temperature for solving the Boltzmann equation
     * @note
     * @param  T_BEGIN: initial temperature
     * @param  T_END: final temperature
     * @retval None
     */
    void Set_T_Range(REAL T_BEGIN, REAL T_END);

    /**
     * @brief  Setting z = scale/T range for solving the Boltzmann equation
     * @note
     * @param  Z_BEGIN: initial z value
     * @param  Z_END: final z value
     * @retval None
     */
    void Set_Z_Range(REAL Z_BEGIN, REAL Z_END);

    /**
     * @brief  Solve the Boltzmann equation
     * @note
     * @param  step_size: initial step size
     * @param  eps_rel: tolerance for max relative error
     * @retval
     */
    RungeKutta::STATUS Solve(REAL step_size, REAL eps_rel = 1e-4);

    /**
     * @brief  Dump current solution of the BE to file
     * @note
     * @param  filename: file name for output
     * @retval None
     */
    void Dump_Solution(std::string filename);
};

}  // namespace EvoEMD

#endif  //_BOLTZMANN_EQUATION_H_
