#include "EvoEMD/BoltzmannEquation.h"

#include <cmath>
#include <fstream>
#include <iomanip>

#include "EvoEMD/EffDOF.h"
#include "EvoEMD/ProcessBase.h"
#include "EvoEMD/spdlog_wrapper.h"

namespace EvoEMD {

Boltzmann_Equation::Boltzmann_Equation(Parameter_Base *ptr_scale_)
    : pf(EvoEMD::Particle_Factory::Get_Particle_Factory()),
      hf(EvoEMD::Hubble_Factory::Get_Hubble_Factory()),
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
    // * sH/gestar(beta_R dY/dx + 3(gestar-beta_R gsstar)Y) = CollisionRate
    // * beta_R indicates the dependence of rho_R on scale factor: rho_R~ a^{-beta_R}

    for (int id = 0; id < DOF; id++) {
        poi_ptrs[id]->Yield = y[id];
        poi_ptrs[id]->Delta_Yield_Ratio = delta_y_ratio[id];
    }

    REAL z = exp(x);
    REAL T = scale / z;
    Hubble_Base *hh = hf.Get_Hubble_Calculator();
    double dlna_dlnT = hh->Get_dlna_dlnT_at_T(T);
    REAL gsstar_at_T = f_gs_star(T);
    REAL gs_at_T = f_gs(T);
    REAL entropy_at_T = 2 * M_PI * M_PI / 45.0 * gs_at_T * pow(T, 3);
    REAL HatT = hh->Get_Hubble_at_T(T);
    REAL res = 0;

    Particle_Base *pp = poi_ptrs[i];
    std::set<Process *> sp = pp->Get_Process();
    for (auto &&proc_ptr : sp) {
        res += proc_ptr->Get_Collision_Rate(T) * proc_ptr->Get_Offset(T, pp->Get_PID());
    }
    res *= -dlna_dlnT / entropy_at_T / HatT;
    res += 3.0 * (gsstar_at_T + dlna_dlnT) * y[i];
    return res;
}

VD Boltzmann_Equation::dYdX(REAL x, const VD &y, const VD &delta_y_ratio) {
    // * x = log(z), z = scale/T
    // * Then the Boltzmann equation for particle Yield evolution is
    // * T^3H(beta_T dY/dx + 3(1-beta_R)Y) = CollisionRate
    // * beta_T indicates the dependence of T on scale factor: T~ a^{-beta_T}

    for (int i = 0; i < DOF; i++) {
        poi_ptrs[i]->Yield = y[i];
        poi_ptrs[i]->Delta_Yield_Ratio = delta_y_ratio[i];
    }

    REAL z = exp(x);
    REAL T = scale / z;
    Hubble_Base *hh = hf.Get_Hubble_Calculator();
    double dlna_dlnT = hh->Get_dlna_dlnT_at_T(T);
    REAL gsstar_at_T = f_gs_star(T);
    REAL gs_at_T = f_gs(T);
    REAL entropy_at_T = 2 * M_PI * M_PI / 45.0 * gs_at_T * pow(T, 3);
    REAL HatT = hh->Get_Hubble_at_T(T);
    // hs->Print();
    VD res(DOF, 0);
    for (int i = 0; i < DOF; i++) {
        Particle_Base *pp = poi_ptrs[i];
        std::set<Process *> sp = pp->Get_Process();
        for (auto &&proc_ptr : sp) {
            REAL cr = proc_ptr->Get_Collision_Rate(T);
            REAL coef = proc_ptr->Get_Offset(T, pp->Get_PID());
            SPDLOG_DEBUG_FILE("COMPONENT-{}, Rate: {:+9.8e}, COEF: {:+9.8e}", i, cr, coef);
            res[i] += cr * coef;
        }
    }
    for (int i = 0; i < DOF; i++) {
        res[i] *= -dlna_dlnT / entropy_at_T / HatT;
        res[i] += 3.0 * (gsstar_at_T + dlna_dlnT) * y[i];
    }
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

VD Boltzmann_Equation::Get_Yield_at_T_End() { return rk.Get_Solution_Y_END(); }

VD Boltzmann_Equation::Get_Omegah2_at_Today() {
    static const REAL rhocoverh2 = 8.098e-47;       // * = rhoc/h^2 in GeV^4;
    static const REAL T0 = 2.348e-13;               // * in GeV;
    static const REAL s_pre_factor = 1.71682909;    // *
    static const REAL h2T03overrhoc = 1.59851e8;    // * = T0^3/rhocoverh2 in GeV^-1;
    static const REAL h2soverrhoc = 2.744375723e8;  // * = s/rhocoverh2 in GeV^-1

    VD Yend = Get_Yield_at_T_End();
    VD Omegah2(DOF);
    for (int i = 0; i < DOF; i++) {
        Omegah2[i] = Yend[i] * poi_ptrs[i]->Get_Mass() * h2soverrhoc;
    }
    return Omegah2;
}

}  // namespace EvoEMD
