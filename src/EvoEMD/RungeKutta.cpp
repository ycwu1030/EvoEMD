#include "EvoEMD/RungeKutta.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "EvoEMD/spdlog_wrapper.h"

using namespace std;

namespace EvoEMD {
void ODE_FUNCS::Set_BOUNDARIES(REAL X_BEGIN, REAL X_END, VD boundary, bool relative_to_equilibrium) {
    this->X_BEGIN = X_BEGIN;
    this->X_END = X_END;
    VD YEQ = Yeq(X_BEGIN);
    if (DOF < 0) {
        std::cout << "You forgot to set DOF first, but we guess you want to use: boundary.size() = " << boundary.size()
                  << std::endl;
        Set_DOF(boundary.size());
    }
    if (boundary.size() != DOF) {
        std::cout << "WRONG dimension in boundary: " << boundary.size() << "!=" << DOF << std::endl;
        return;
    }
    if (!relative_to_equilibrium) {
        for (int i = 0; i < DOF; i++) {
            Y_BEGIN[i] = boundary[i];
            Delta_Y_Ratio_BEGIN[i] = 1.0 - Y_BEGIN[i] / YEQ[i];
        }
    } else {
        for (int i = 0; i < DOF; i++) {
            Y_BEGIN[i] = boundary[i] * YEQ[i];
            Delta_Y_Ratio_BEGIN[i] = 1.0 - boundary[i];
        }
    }
}

RungeKutta::RungeKutta()
    : DOF(0),
      derivs(nullptr),
      SOLVED(false),
      SAFETY(0.9),
      POW_GROW(-0.2),
      POW_SHRINK(-0.25),
      TINY_REL_THRESHOLD(1e-60),
      MAXSTEPS(10000) {}

RungeKutta::RungeKutta(ODE_FUNCS *ode)
    : SOLVED(false), SAFETY(0.9), POW_GROW(-0.2), POW_SHRINK(-0.25), TINY_REL_THRESHOLD(1e-60), MAXSTEPS(10000) {
    Set_ODE(ode);
}

void RungeKutta::Set_ODE(ODE_FUNCS *ode) {
    derivs = ode;
    SOLVED = false;
}

void RungeKutta::Clear() {
    SPDLOG_INFO_FILE("Clear X, Y and dYdX in RungeKutta Solver.");
    _X.clear();
    _Y.clear();
    _Delta_Y_Ratio.clear();
    _Yeq.clear();
    _dYdX.clear();
}

void RungeKutta::INIT() {
    SPDLOG_INFO_FILE("Initialize the RungeKutta Solver.");
    // * Clear whatever we have in our system
    Clear();

    // * Getting information from the ODE system
    DOF = derivs->Get_DOF();
    X_BEGIN = derivs->Get_X_BEGIN();
    X_END = derivs->Get_X_END();
    Y_BEGIN = derivs->Get_Y_BEGIN();
    Delta_Y_Ratio_BEGIN = derivs->Get_Delta_Y_Ratio_BEGIN();
    ZERO_THRESHOLD = Y_BEGIN * TINY_REL_THRESHOLD;

    // * Initialize the system at X_BEGIN:
    _X.push_back(X_BEGIN);
    _Yeq.push_back(derivs->Yeq(X_BEGIN));
    _Y.push_back(Y_BEGIN);
    _Delta_Y_Ratio.push_back(Delta_Y_Ratio_BEGIN);

    _dYdX.push_back((*derivs)(_X[0], _Y[0], _Delta_Y_Ratio[0]));

    SPDLOG_DEBUG_FILE("Initialize the RungeKutta with:");
    SPDLOG_DEBUG_FILE("DOF: {}.", DOF);
    SPDLOG_DEBUG_FILE("ODE Range: [{:+9.8e}, {:+9.8e}].", X_BEGIN, X_END);
    SPDLOG_DEBUG_FILE("Boundary Condition at x = {:+9.8e}:", X_BEGIN);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e}, ZERO_THRESHOLD[{0}] = {3:+9.8e};", i, _Y[0][i],
                          _dYdX[0][i], ZERO_THRESHOLD[i]);
    }
}

bool RungeKutta::RK4_SingleStep(const REAL x_cur, const VD &y_cur, const VD &dydx_cur, const VD &delta_y_ratio_cur,
                                const REAL step_size, VD &y_next, VD &delta_y_ratio_next) {
    SPDLOG_INFO_FILE("4th Order RK with stepsize = {:+9.8e}.", step_size);
    SPDLOG_INFO_FILE("Starting at: ");
    SPDLOG_INFO_FILE("  X = {:+9.8e}.", x_cur);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e};", i, y_cur[i], dydx_cur[i]);
    }
    REAL half_step = step_size / 2.0;

    // * 1. Using dx*dy/dx
    VD dY_Step1 = step_size * dydx_cur;
    SPDLOG_DEBUG_FILE("  4th Order RK Step-1:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY1[{}] = {:+9.8e};", i, dY_Step1[i]);
    }

    // * 2. Using dx/2 and dY_Step1/2
    VD Yeq_half_step = derivs->Yeq(x_cur + half_step);
    VD delta_y_ratio_step2 = 1.0 - (y_cur + dY_Step1 / 2.0) / Yeq_half_step;
    VD dY_Step2 = step_size * (*derivs)(x_cur + half_step, y_cur + dY_Step1 / 2.0, delta_y_ratio_step2);
    SPDLOG_DEBUG_FILE("  4th Order RK Step-2:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY2[{}] = {:+9.8e};", i, dY_Step2[i]);
    }

    // * 3. Using dx/2 and dY_Step2/2
    VD delta_y_ratio_step3 = 1.0 - (y_cur + dY_Step2 / 2.0) / Yeq_half_step;
    VD dY_Step3 = step_size * (*derivs)(x_cur + half_step, y_cur + dY_Step2 / 2.0, delta_y_ratio_step3);
    SPDLOG_DEBUG_FILE("  4th Order RK Step-3:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY3[{}] = {:+9.8e};", i, dY_Step3[i]);
    }

    // * 4. Using dx and dY_Step3
    VD Yeq_full_step = derivs->Yeq(x_cur + step_size);
    VD delta_y_ratio_step4 = 1.0 - (y_cur + dY_Step3) / Yeq_full_step;
    VD dY_Step4 = step_size * (*derivs)(x_cur + step_size, y_cur + dY_Step3, delta_y_ratio_step4);
    SPDLOG_DEBUG_FILE("  4th Order RK Step-4:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY4[{}] = {:+9.8e};", i, dY_Step4[i]);
    }

    // * Combine above 4 steps:
    y_next = y_cur + dY_Step1 / 6.0 + dY_Step2 / 3.0 + dY_Step3 / 3.0 + dY_Step4 / 6.0;
    SPDLOG_INFO_FILE("4th Order RK with stepsize = {:+9.8e}.", step_size);
    SPDLOG_INFO_FILE("Ending at: ");
    SPDLOG_INFO_FILE("  X = {:+9.8e}.", x_cur + step_size);
    std::vector<bool> status = derivs->Is_Thermalized();
    bool comp = false;
    for (int i = 0; i < DOF; i++) {
        SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e};", i, y_next[i]);
        if (status[i] && y_next[i] < Yeq_full_step[i]) {
            if (y_next[i] < Yeq_full_step[i]) {
                SPDLOG_INFO_FILE("    From Y[{0}] = Yeq[{0}] = {1:+9.8e} to Y[{0}] = {2:+9.8e}", i, y_cur[i],
                                 y_next[i]);
                SPDLOG_INFO_FILE("    Y[{0}] = {1:+9.8e} < Yeq[{0}] = {2:+9.8e}", i, y_next[i], Yeq_full_step[i]);
                SPDLOG_INFO_FILE("    Compensate back to Yeq");
                comp = true;
                y_next[i] = Yeq_full_step[i];
                status[i] = true;
                delta_y_ratio_next[i] = 0;
            } else {
                status[i] = false;
                delta_y_ratio_next[i] = 1.0 - y_next[i] / Yeq_full_step[i];
            }
        } else {
            // * Checking the situation at first order
            if (delta_y_ratio_cur[i] > 0) {
                SPDLOG_INFO_FILE("    Start from Y[{0}] = {1:+9.8e} < Yeq[{0}]", i, y_cur[i]);
                if (y_cur[i] + dY_Step1[i] > Yeq_full_step[i]) {
                    SPDLOG_INFO_FILE("    Changing too fast, compensate to Yeq.");
                    comp = true;
                    y_next[i] = Yeq_full_step[i];
                    status[i] = true;
                    delta_y_ratio_next[i] = 0;
                } else {
                    status[i] = false;
                    delta_y_ratio_next[i] = 1.0 - y_next[i] / Yeq_full_step[i];
                }
            } else {
                SPDLOG_INFO_FILE("    Start from Y[{0}] = {1:+9.8e} > Yeq[{0}]", i, y_cur[i]);
                if (y_cur[i] + dY_Step1[i] < Yeq_full_step[i]) {
                    SPDLOG_INFO_FILE("    Changing too fast, compensate to Yeq.");
                    comp = true;
                    y_next[i] = Yeq_full_step[i];
                    status[i] = true;
                    delta_y_ratio_next[i] = 0;
                } else {
                    status[i] = false;
                    delta_y_ratio_next[i] = 1.0 - y_next[i] / Yeq_full_step[i];
                }
            }
        }
    }
    derivs->Update_Thermal_Status(status);
    return comp;
}

bool RungeKutta::RKQC_SingleStep(REAL &x, VD &y, VD &yeq, VD &dydx, VD &delta_y_ratio, const REAL step_size_guess,
                                 const REAL eps, const VD &Y_Scale, REAL &step_size_did, REAL &step_size_further) {
    SPDLOG_INFO_FILE("RKQC: ");
    SPDLOG_INFO_FILE("Starting with: ");
    SPDLOG_INFO_FILE("Guess step size: dx = {:+9.8e},", step_size_guess);
    SPDLOG_INFO_FILE("  X = {:+9.8e}.", x);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e}, Yeq[{0}] = {3:+9.8e};", i, y[i], dydx[i],
                         yeq[i]);
    }
    SPDLOG_INFO_FILE("  Guessed Step size dx = {:+9.8e}", step_size_guess);
    // * Cache the initial points
    REAL x_cache = x;
    VD y_cache = y;
    VD yeq_cache = yeq;
    VD dy_cache = dydx;
    VD delta_y_ratio_cache = delta_y_ratio;
    REAL step_size = step_size_guess;
    REAL step_size_temp;
    REAL half_step_size;

    VD y_temp(DOF, 0);
    VD delta_y_ratio_temp(DOF, 0);
    VD Delta_y(DOF, 0);
    VD error_temp(DOF, 0);
    std::vector<bool> status = derivs->Is_Thermalized();

    REAL error_max = 0;
    REAL min_step_size = 1e-5 * step_size_guess;
    SPDLOG_DEBUG_FILE("Minimal allowed step size is {:+9.8e};", min_step_size);
    REAL max_step_size;
    int trials = 0;
    bool comp_two_step;
    bool comp_single_step;
    bool comp;
    while (true) {
        ++trials;
        SPDLOG_DEBUG_FILE("{}-th Trial for RKQC with stepsize = {:+9.8e}", trials, step_size);
        // * Take two half steps
        half_step_size = step_size / 2.0;
        derivs->Update_Thermal_Status(status);
        RK4_SingleStep(x_cache, y_cache, dy_cache, delta_y_ratio_cache, half_step_size, y_temp, delta_y_ratio_temp);
        x = x_cache + half_step_size;
        dydx = (*derivs)(x, y_temp, delta_y_ratio_temp);
        comp_two_step = RK4_SingleStep(x, y_temp, dydx, delta_y_ratio_temp, half_step_size, y, delta_y_ratio);
        SPDLOG_DEBUG_FILE("Take two half steps, reaching:");
        SPDLOG_DEBUG_FILE("  X = {:+9.8e};", x + half_step_size);
        for (int i = 0; i < DOF; i++) {
            if (status[i] && y[i] < ZERO_THRESHOLD[i]) {
                SPDLOG_WARN_FILE("  REACH ZERO THRESHOLD AS Y[{0}] = {1:+9.8e} < {2:+9.8e};", i, y[i],
                                 ZERO_THRESHOLD[i]);
                y[i] = 0;
            }
            SPDLOG_DEBUG_FILE("  Y[{}] = {:+9.8e};", i, y[i]);
        }

        // * Take one full step
        derivs->Update_Thermal_Status(status);
        comp_single_step =
            RK4_SingleStep(x_cache, y_cache, dy_cache, delta_y_ratio_cache, step_size, y_temp, delta_y_ratio_temp);
        bool drop_slow = true;
        for (int i = 0; i < DOF; i++) {
            REAL slop = std::fabs(y_cache[i] / y_temp[i]);
            drop_slow *= (slop < 5);
        }
        x = x_cache + step_size;
        SPDLOG_DEBUG_FILE("Take one full step, reaching:");
        SPDLOG_DEBUG_FILE("  X = {:+9.8e};", x);
        for (int i = 0; i < DOF; i++) {
            if (status[i] && y_temp[i] < ZERO_THRESHOLD[i]) {
                SPDLOG_WARN_FILE("  REACH ZERO THRESHOLD AS Y[{0}] = {1:+9.8e} < {2:+9.8e};", i, y_temp[i],
                                 ZERO_THRESHOLD[i]);
                y_temp[i] = 0;
            }
            SPDLOG_DEBUG_FILE("  Y[{}] = {:+9.8e};", i, y_temp[i]);
        }
        comp = comp_single_step && comp_two_step;
        // * Check the difference between above two methods

        if (comp && !drop_slow) {
            // * We compensate to the Equilibrium, and the yeq drops too fast
            // * In order not to miss the freeze-out point, we have to shrink the stepsize;
            step_size_temp = step_size;
            step_size = step_size_temp / 2.0;
            SPDLOG_WARN_FILE(
                "  We compensate to Equilibrium Value, but it drops too fast, we need to shrink the step size to be "
                "safe:  {:+9.8e} -> {:+9.8e}",
                step_size_temp, step_size);
            continue;
        }
        Delta_y = y - y_temp;
        error_temp = fabs(Delta_y / Y_Scale);
        error_max = *max_element(error_temp.begin(), error_temp.end());
        error_max /= eps;
        for (int i = 0; i < DOF; i++) {
            SPDLOG_DEBUG_FILE("Calc Error: DeltaY[{0}] = {1:+9.8e},  YScale[{0}] = {2:+9.8e}, error[{0}] = {3:+9.8e}",
                              i, Delta_y[i], Y_Scale[i], error_temp[i]);
        }

        SPDLOG_DEBUG_FILE("Error max: {:+9.8e};", error_max);
        if (!comp_two_step && comp_single_step) {
            step_size_temp = step_size;
            step_size = step_size_temp / 2.0;
            SPDLOG_INFO_FILE(
                "We are crossing the critical point, need to shrink the stepsize in order to probe it: {:+9.8e} -> "
                "{:+9.8e}",
                step_size_temp, step_size);
            continue;
        }
        if (error_max <= 1.0) {
            // * The error is acceptable, we will proceed, and even try to enlarge the step size
            step_size_did = step_size;
            step_size_temp = SAFETY * step_size * exp(POW_GROW * log(error_max + 1e-5));
            max_step_size = 4 * step_size_did;
            step_size_further = (comp_single_step) ? 1.5 * step_size_did : min(step_size_temp, max_step_size);
            SPDLOG_DEBUG_FILE("Error is small, enlarge the step size to: {:+9.8e};", step_size_further);
            break;
        }
        // * Otherwise, the error is larger than our tolarance, we need to shrink the step size to improve the error
        // * But, we need to first check whether the error is tooooo large. If so restrict it to a reasonable value
        // * Actually 1e10 is already not that reasonable, but at least, it should enough to avoid too small step_size
        error_max = error_max < 1e10 ? error_max : 1e10;
        step_size_temp = SAFETY * step_size * exp(POW_SHRINK * log(error_max));
        SPDLOG_DEBUG_FILE("Error is large, try to shrink the step size to: {:+9.8e};", step_size_temp);
        if (std::fabs(step_size_temp) < std::fabs(min_step_size)) {
            // * If after shrink, the step size is too small, then we don't shrink it any further, and accept current
            // possible large error.
            step_size_did = step_size;
            step_size_further = step_size;
            SPDLOG_WARN_FILE("Might reach the smallest step size");
            break;
        }
        step_size = step_size_temp;  // * Change the step size, do the RK4 again.
    }
    y = y + Delta_y / 15;
    yeq = derivs->Yeq(x);
    if (comp) {
        for (int i = 0; i < DOF; i++) {
            delta_y_ratio[i] = 0;
        }
    } else {
        for (int i = 0; i < DOF; i++) {
            delta_y_ratio[i] = 1.0 - y[i] / yeq[i];
        }
    }
    dydx = (*derivs)(x, y, delta_y_ratio);
    SPDLOG_DEBUG_FILE("RKQC reaching: ");
    SPDLOG_DEBUG_FILE("  X = {:+9.8e};", x);
    SPDLOG_DEBUG_FILE("  dx = {:+9.8e};", step_size_did);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE(" Y[{0}] = {1:+9.8e}, dYdX = {2:+9.8e};", i, y[i], dydx[i]);
    }
    return comp;
}

RungeKutta::STATUS RungeKutta::Solve(REAL step_start, REAL eps) {
    // * Initialize all relevant quantities
    SPDLOG_INFO_FILE("Start to solving the ODE with RungeKutta Method.");
    INIT();

    // * The initial point
    double x = _X[0];
    VD y = _Y[0];
    VD yeq = _Yeq[0];
    VD dydx = _dYdX[0];
    VD delta_y_ratio = _Delta_Y_Ratio[0];
    VD Y_Scale(DOF);

    double step_size = (X_END > X_BEGIN) ? std::fabs(step_start) : -std::fabs(step_start);
    double step_size_did;
    double step_size_next;

    for (int nstp = 0; nstp < MAXSTEPS; nstp++) {
        SPDLOG_INFO_FILE("Step-{}:", nstp + 1);
        SPDLOG_INFO_FILE("  X = {:+9.8e};", x);
        Y_Scale = fabs(y) + fabs(dydx * step_size);
        for (int i = 0; i < DOF; i++) {
            // Y_Scale[i] = min(1.0, Y_Scale[i]);
            Y_Scale[i] = max(1e-50, Y_Scale[i]);
            SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e}, YScale[{0}] = {3:+9.8e};", i, y[i], dydx[i],
                             Y_Scale[i]);
        }

        // * Check whether the step size is too large that we already pass the end point
        if ((step_size > 0 && x + step_size > X_END) || (step_size < 0 && x + step_size < X_END)) {
            step_size = X_END - x;
        }

        // * One step forward
        bool comp =
            RKQC_SingleStep(x, y, yeq, dydx, delta_y_ratio, step_size, eps, Y_Scale, step_size_did, step_size_next);

        _X.push_back(x);
        _Y.push_back(y);
        _Yeq.push_back(yeq);
        _dYdX.push_back(dydx);
        _Delta_Y_Ratio.push_back(delta_y_ratio);

        // * Check if we reach the end point
        if ((x - X_END) * (X_END - X_BEGIN) >= 0) {
            SOLVED = true;
            SPDLOG_INFO_FILE("Solve the ODE successfully.");
            return SUCCESS;
        }

        // * If not reaching the end, we adapt the stepsize, but need to control it either not too large or too small.
        step_size = step_size_next;
        if (std::fabs(step_size_next) > 1e-1 * std::fabs(X_END - X_BEGIN)) {
            step_size = 1e-1 * (X_END - X_BEGIN);
        }
        if (std::fabs(step_size_next) < 1e-5 * std::fabs(X_END - X_BEGIN)) {
            step_size = 1e-5 * (X_END - X_BEGIN);
        }
    }
    SOLVED = true;
    SPDLOG_WARN_FILE("Solve the ODE with too many steps.");
    return TOOMANYSTEP;
}

void RungeKutta::Dump_Solution(string filename) {
    ofstream output(filename.c_str());
    output << "x\t";
    for (size_t i = 0; i < DOF; i++) {
        output << "y_" << i << "\t"
               << "yeq_" << i << "\t";
    }
    for (size_t i = 0; i < DOF; i++) {
        output << "dydx_" << i << "\t";
    }
    output << endl;
    output << scientific;
    output << showpos;
    output << setprecision(8);
    for (size_t i = 0; i < _X.size(); i++) {
        output << _X[i] << "\t";
        for (size_t j = 0; j < DOF; j++) {
            output << _Y[i][j] << "\t" << _Yeq[i][j] << "\t";
        }
        for (size_t j = 0; j < DOF; j++) {
            output << _dYdX[i][j] << "\t";
        }
        output << endl;
    }
    output.close();
}
}  // namespace EvoEMD
