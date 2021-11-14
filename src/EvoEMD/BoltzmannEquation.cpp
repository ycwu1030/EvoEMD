#include "EvoEMD/BoltzmannEquation.h"

#include <cmath>
#include <fstream>
#include <iomanip>

#include "EvoEMD/ProcessBase.h"
#include "EvoEMD/spdlog_wrapper.h"

namespace EvoEMD {

Boltzmann_Equation::Boltzmann_Equation(Parameter_Base *ptr_scale_)
    : pf(EvoEMD::Particle_Factory::Get_Particle_Factory()),
      hh(EvoEMD::Hubble_History::Get_Hubble_History()),
      ptr_scale(ptr_scale_),
      rk(this) {
    std::set<int> POI = pf.Get_POI();
    Set_DOF(POI.size());
    std::set<int>::iterator iter = POI.begin();
    for (; iter != POI.end(); iter++) {
        poi_pids.push_back(*iter);
        Particle_Base *pp = pf.Get_Particle(*iter);
        poi_ptrs.push_back(pp);
        poi_names.push_back(pp->Get_Name());
    }
    Setup_Scale();
}

void Boltzmann_Equation::Setup_Scale() {
    if (!ptr_scale) {
        scale = 1000.0;
    } else {
        scale = ptr_scale->Get_Value();
    }
}

VB Boltzmann_Equation::Should_be_Thermalized(REAL x, const VD &y, const VD &delta_y_ratio) {
    VB res(DOF, false);
    VD yeq = Yeq(x);
    for (int i = 0; i < DOF; i++) {
        Particle_Base *pp = poi_ptrs[i];
        bool pseudo = pp->Is_Pseudo();
        if (!pseudo && fabs(delta_y_ratio[i]) < 1e-3) {
            // * If this component is a basic particle, and it is close to its equilibrium number density
            // * We use 1e-4 step size to test whether it is possible to thermalized it
            // * when deviated from equilibrium by 1e-3
            VD y_test = y;
            VD delta_y_ratio_test = delta_y_ratio;
            y_test[i] = (1.0 + 1e-3) * yeq[i];
            delta_y_ratio_test[i] = -1e-3;
            REAL dydx_up = dYidX(i, x, y_test, delta_y_ratio_test);
            y_test[i] = (1.0 - 1e-3) * yeq[i];
            delta_y_ratio_test[i] = 1e-3;
            REAL dydx_down = dYidX(i, x, y_test, delta_y_ratio_test);
            SPDLOG_DEBUG_FILE("Component-{}, dydx_up = {:+9.8e}, dydx_down = {:+9.8e}, yeq = {:+9.8e}", i, dydx_up,
                              dydx_down, yeq[i]);
            if (dydx_up < 0 && fabs(dydx_up) * 1e-4 > 1e-3 * yeq[i] && dydx_down > 0 &&
                dydx_down * 1e-4 > 1e-3 * yeq[i]) {
                res[i] = true;
            }
        }

        if (fabs(delta_y_ratio[i]) > 1e-3) {
            REAL dydx = dYidX(i, x, y, delta_y_ratio);
            if ((y[i] + dydx * 1e-4 - yeq[i]) * (y[i] - yeq[i]) < 0) {
                res[i] = true;
            }
        }
    }
    return res;
}

REAL Boltzmann_Equation::dYidX(int i, REAL x, const VD &y, const VD &delta_y_ratio) {
    // * x = log(z), z = scale/T
    // * Then the Boltzmann equation for particle Yield evolution is
    // * T^3H(beta_T dY/dx + 3(1-beta_T)Y) = CollisionRate
    // * beta_T indicates the dependence of T on scale factor: T~ a^{-beta_T}

    for (int id = 0; id < DOF; id++) {
        poi_ptrs[id]->Yield = y[id];
        poi_ptrs[id]->Delta_Yield_Ratio = delta_y_ratio[id];
    }

    REAL z = exp(x);
    REAL T = scale / z;
    int hubble_period_id = hh.Get_Period_ID_at_T(T);
    Hubble_For_Single_Period *hs = hh[hubble_period_id];
    double beta_T = hs->Get_beta_T();
    REAL HatT = hs->Get_Hubble_at_T(T);
    bool isentropic = hs->Is_Isentropic();
    REAL res = 0;

    Particle_Base *pp = poi_ptrs[i];
    std::set<Process *> sp = pp->Get_Process();
    for (auto &&proc_ptr : sp) {
        // std::cout << "Collision Rate @ T=" << T << " is " << proc_ptr->Get_Collision_Rate(T) << std::endl;
        // std::cout << "Coefficient @ T=" << T << " is " << proc_ptr->Get_Yield_Coeff(T, pp->Get_PID()) <<
        // std::endl;
        res += proc_ptr->Get_Collision_Rate(T) * proc_ptr->Get_Yield_Coeff(T, pp->Get_PID());
    }
    res *= 1.0 / pow(T, 3) / HatT;
    if (!isentropic) {
        // that beta_T != 1
        res -= 3.0 * (1.0 - beta_T) * y[i];
        res /= beta_T;
    }  // else beta_T == 1, no further action needed
    return res;
}

VD Boltzmann_Equation::dYdX(REAL x, const VD &y, const VD &delta_y_ratio) {
    // * x = log(z), z = scale/T
    // * Then the Boltzmann equation for particle Yield evolution is
    // * T^3H(beta_T dY/dx + 3(1-beta_T)Y) = CollisionRate
    // * beta_T indicates the dependence of T on scale factor: T~ a^{-beta_T}

    for (int i = 0; i < DOF; i++) {
        poi_ptrs[i]->Yield = y[i];
        poi_ptrs[i]->Delta_Yield_Ratio = delta_y_ratio[i];
    }

    REAL z = exp(x);
    REAL T = scale / z;
    int hubble_period_id = hh.Get_Period_ID_at_T(T);
    Hubble_For_Single_Period *hs = hh[hubble_period_id];
    double beta_T = hs->Get_beta_T();
    REAL HatT = hs->Get_Hubble_at_T(T);
    bool isentropic = hs->Is_Isentropic();
    // hs->Print();
    VD res(DOF, 0);
    for (int i = 0; i < DOF; i++) {
        Particle_Base *pp = poi_ptrs[i];
        std::set<Process *> sp = pp->Get_Process();
        for (auto &&proc_ptr : sp) {
            // std::cout << "Collision Rate @ T=" << T << " is " << proc_ptr->Get_Collision_Rate(T) << std::endl;
            // std::cout << "Coefficient @ T=" << T << " is " << proc_ptr->Get_Yield_Coeff(T, pp->Get_PID()) <<
            // std::endl;
            REAL cr = proc_ptr->Get_Collision_Rate(T);
            REAL coef = proc_ptr->Get_Yield_Coeff(T, pp->Get_PID());
            SPDLOG_DEBUG_FILE("COMPONENT-{}, Rate: {:+9.8e}, COEF: {:+9.8e}", i, cr, coef);
            res[i] += cr * coef;
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
    // std::cout << "RES = " << res << std::endl;
    return res;
}

VD Boltzmann_Equation::Yeq(REAL x) {
    REAL z = exp(x);
    REAL T = scale / z;
    VD res(DOF, 0);
    for (int i = 0; i < DOF; i++) {
        Particle_Base *pp = poi_ptrs[i];
        REAL yeq = pp->Get_Equilibrium_Yield_at_T(T);
        res[i] = yeq;
    }
    return res;
}

VB Boltzmann_Equation::Is_Thermalized() {
    VB res(DOF);
    for (int i = 0; i < DOF; i++) {
        res[i] = poi_ptrs[i]->Get_Init_Thermal_Status();
    }
    return res;
}

VB Boltzmann_Equation::Can_be_Negative() {
    VB res(DOF);
    for (int i = 0; i < DOF; i++) {
        // * If it is pseudo, then the number density/Yield can be negative.
        res[i] = poi_ptrs[i]->Is_Pseudo();
    }
    return res;
}

void Boltzmann_Equation::Set_X_Range(REAL X_BEGIN, REAL X_END) {
    this->X_BEGIN = X_BEGIN;
    this->X_END = X_END;
    Y_BEGIN = Yeq(X_BEGIN);
    std::vector<bool> status = Is_Thermalized();
    for (int i = 0; i < DOF; i++) {
        if (status[i]) {
            Delta_Y_Ratio_BEGIN[i] = 0.0;
        } else {
            Y_BEGIN[i] = 0;
            Delta_Y_Ratio_BEGIN[i] = 1.0;
        }
    }
}

void Boltzmann_Equation::Set_T_Range(REAL T_BEGIN, REAL T_END) {
    this->T_BEGIN = T_BEGIN;
    this->T_END = T_END;
}

void Boltzmann_Equation::Set_Z_Range(REAL Z_BEGIN, REAL Z_END) {
    Setup_Scale();
    SPDLOG_INFO_FILE("Setting Z range with scale = {:+9.8e}", scale);
    T_BEGIN = scale / Z_BEGIN;
    T_END = scale / Z_END;
}

RungeKutta::STATUS Boltzmann_Equation::Solve(REAL step_size, REAL eps_rel) {
    Setup_Scale();
    REAL x_begin = log(scale / T_BEGIN);
    REAL x_end = log(scale / T_END);
    Set_X_Range(x_begin, x_end);
    RungeKutta::STATUS st = rk.Solve(step_size, eps_rel);
    Process_Factory::Clean_Cache();
    return st;
}

void Boltzmann_Equation::Dump_Solution(std::string filename) {
    std::ofstream output(filename.c_str());
    VD X = rk.Get_Solution_X();
    VVD Y = rk.Get_Solution_Y();
    VVD dYdX = rk.Get_Solution_dYdX();
    VVD Yeq = rk.Get_Solution_Yeq();
    output << "z\t";
    for (size_t i = 0; i < DOF; i++) {
        output << "y_" << i << "\t";
    }
    for (size_t i = 0; i < DOF; i++) {
        output << "yeq_" << i << "\t";
    }
    for (size_t i = 0; i < DOF; i++) {
        output << "dydz_" << i << "\t";
    }
    output << "scale" << std::endl;
    output << std::scientific;
    output << std::showpos;
    output << std::setprecision(8);
    for (size_t i = 0; i < X.size(); i++) {
        REAL z = exp(X[i]);
        output << z << "\t";
        for (size_t j = 0; j < DOF; j++) {
            output << Y[i][j] << "\t";
        }
        for (size_t j = 0; j < DOF; j++) {
            output << Yeq[i][j] << "\t";
        }
        for (size_t j = 0; j < DOF; j++) {
            output << dYdX[i][j] / z << "\t";  // * from dY/dX -> dY/dZ
        }
        output << scale << std::endl;
    }
    output.close();
}

}  // namespace EvoEMD
