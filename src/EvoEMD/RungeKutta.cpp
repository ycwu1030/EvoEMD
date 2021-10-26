#include "EvoEMD/RungeKutta.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "EvoEMD/spdlog_wrapper.h"

using namespace std;

namespace EvoEMD {
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
    BOUNDARY_AT_BEGIN = derivs->Get_BOUNDARY_CONDITION();
    ZERO_THRESHOLD = BOUNDARY_AT_BEGIN * TINY_REL_THRESHOLD;

    // * Initialize the system at X_BEGIN:
    _X.push_back(X_BEGIN);
    _Y.push_back(BOUNDARY_AT_BEGIN);
    _Yeq.push_back(derivs->Yeq(X_BEGIN));
    _dYdX.push_back((*derivs)(X_BEGIN, BOUNDARY_AT_BEGIN));

    SPDLOG_DEBUG_FILE("Initialize the RungeKutta with:");
    SPDLOG_DEBUG_FILE("DOF: {}.", DOF);
    SPDLOG_DEBUG_FILE("ODE Range: [{:+9.8e}, {:+9.8e}].", X_BEGIN, X_END);
    SPDLOG_DEBUG_FILE("Boundary Condition at x = {:+9.8e}:", X_BEGIN);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e}, ZERO_THRESHOLD[{0}] = {3:+9.8e};", i,
                          BOUNDARY_AT_BEGIN[i], _dYdX[0][i], ZERO_THRESHOLD[i]);
    }
}

bool RungeKutta::RK4_SingleStep(const REAL x_cur, const VD &y_cur, const VD &dy_cur, const REAL step_size, VD &y_next) {
    SPDLOG_INFO_FILE("4th Order RK with stepsize = {:+9.8e}.", step_size);
    SPDLOG_INFO_FILE("Starting at: ");
    SPDLOG_INFO_FILE("  X = {:+9.8e}.", x_cur);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e};", i, y_cur[i], dy_cur[i]);
    }
    REAL half_step = step_size / 2.0;

    // * 1. Using dx*dy/dx
    VD dY_Step1 = step_size * dy_cur;
    SPDLOG_DEBUG_FILE("  4th Order RK Step-1:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY1[{}] = {:+9.8e};", i, dY_Step1[i]);
    }

    // * 2. Using dx/2 and dY_Step1/2
    VD dY_Step2 = step_size * (*derivs)(x_cur + half_step, y_cur + dY_Step1 / 2.0);
    SPDLOG_DEBUG_FILE("  4th Order RK Step-2:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY2[{}] = {:+9.8e};", i, dY_Step2[i]);
    }

    // * 3. Using dx/2 and dY_Step2/2
    VD dY_Step3 = step_size * (*derivs)(x_cur + half_step, y_cur + dY_Step2 / 2.0);
    SPDLOG_DEBUG_FILE("  4th Order RK Step-3:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY3[{}] = {:+9.8e};", i, dY_Step3[i]);
    }

    // * 4. Using dx and dY_Step3
    VD dY_Step4 = step_size * (*derivs)(x_cur + step_size, y_cur + dY_Step3);
    SPDLOG_DEBUG_FILE("  4th Order RK Step-4:");
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE("    dY4[{}] = {:+9.8e};", i, dY_Step4[i]);
    }

    // * Combine above 4 steps:
    y_next = y_cur + dY_Step1 / 6.0 + dY_Step2 / 3.0 + dY_Step3 / 3.0 + dY_Step4 / 6.0;
    SPDLOG_INFO_FILE("4th Order RK with stepsize = {:+9.8e}.", step_size);
    SPDLOG_INFO_FILE("Ending at: ");
    SPDLOG_INFO_FILE("  X = {:+9.8e}.", x_cur + step_size);
    VD Yeq = derivs->Yeq(x_cur + step_size);
    bool comp = false;
    for (int i = 0; i < DOF; i++) {
        SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e};", i, y_next[i]);
        if (y_next[i] < Yeq[i]) {
            SPDLOG_INFO_FILE("    Y[{0}] = {1:+9.8e} < Yeq[{0}] = {2:+9.8e}", i, y_next[i], Yeq[i]);
            SPDLOG_INFO_FILE("    Compensate back to Yeq");
            comp = true;
            y_next[i] = Yeq[i];
        }
    }
    return comp;
}

void RungeKutta::RKQC_SingleStep(REAL &x, VD &y, VD &dy, const REAL step_size_guess, const REAL eps, const VD &Y_Scale,
                                 REAL &step_size_did, REAL &step_size_further) {
    SPDLOG_INFO_FILE("RKQC: ");
    SPDLOG_INFO_FILE("Starting with: ");
    SPDLOG_INFO_FILE("Guess step size: dx = {:+9.8e},", step_size_guess);
    SPDLOG_INFO_FILE("  X = {:+9.8e}.", x);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e};", i, y[i], dy[i]);
    }
    // * Cache the initial points
    REAL x_cache = x;
    VD y_cache = y;
    VD dy_cache = dy;
    REAL step_size = step_size_guess;
    REAL step_size_temp;
    REAL half_step_size;

    VD y_temp;
    VD Delta_y;
    VD error_temp;

    REAL error_max = 0;
    REAL min_step_size = 1e-5 * step_size_guess;
    SPDLOG_DEBUG_FILE("Minimal allowed step size is {:+9.8e};", min_step_size);
    REAL max_step_size;
    int trials = 0;
    while (true) {
        bool comp = false;
        ++trials;
        SPDLOG_DEBUG_FILE("{}-th Trial for RKQC with stepsize = {:+9.8e}", trials, step_size);
        // * Take two half steps
        half_step_size = step_size / 2.0;
        bool did_comp = RK4_SingleStep(x_cache, y_cache, dy_cache, half_step_size, y_temp);
        comp = comp || did_comp;
        x = x_cache + half_step_size;
        dy = (*derivs)(x, y_temp);
        did_comp = RK4_SingleStep(x, y_temp, dy, half_step_size, y);
        comp = comp || did_comp;
        SPDLOG_DEBUG_FILE("Take two half steps, reaching:");
        SPDLOG_DEBUG_FILE("  X = {:+9.8e};", x + half_step_size);
        for (int i = 0; i < DOF; i++) {
            if (y[i] < ZERO_THRESHOLD[i]) {
                SPDLOG_WARN_FILE("  REACH ZERO THRESHOLD AS Y[{0}] = {1:+9.8e} < {2:+9.8e};", i, y[i],
                                 ZERO_THRESHOLD[i]);
                y[i] = 0;
            }
            SPDLOG_DEBUG_FILE("  Y[{}] = {:+9.8e};", i, y[i]);
        }

        // * Take one full step
        did_comp = RK4_SingleStep(x_cache, y_cache, dy_cache, step_size, y_temp);
        comp = comp || did_comp;
        x = x_cache + step_size;
        SPDLOG_DEBUG_FILE("Take one full step, reaching:");
        SPDLOG_DEBUG_FILE("  X = {:+9.8e};", x);
        for (int i = 0; i < DOF; i++) {
            if (y_temp[i] < ZERO_THRESHOLD[i]) {
                SPDLOG_WARN_FILE("  REACH ZERO THRESHOLD AS Y[{0}] = {1:+9.8e} < {2:+9.8e};", i, y_temp[i],
                                 ZERO_THRESHOLD[i]);
                y_temp[i] = 0;
            }
            SPDLOG_DEBUG_FILE("  Y[{}] = {:+9.8e};", i, y_temp[i]);
        }

        // * Check the difference between above two methods

        Delta_y = y - y_temp;
        error_temp = fabs(Delta_y / Y_Scale);
        error_max = *max_element(error_temp.begin(), error_temp.end());
        error_max /= eps;
        for (int i = 0; i < DOF; i++) {
            SPDLOG_DEBUG_FILE("Calc Error: DeltaY[{0}] = {1:+9.8e},  YScale[{0}] = {2:+9.8e}, error[{0}] = {3:+9.8e}",
                              i, Delta_y[i], Y_Scale[i], error_temp[i]);
        }

        SPDLOG_DEBUG_FILE("Error max: {:+9.8e};", error_max);
        if (error_max <= 1.0) {
            // * The error is acceptable, we will proceed, and even try to enlarge the step size
            step_size_did = step_size;
            step_size_temp = SAFETY * step_size * exp(POW_GROW * log(error_max + 1e-5));
            max_step_size = 4 * step_size_did;
            step_size_further = comp ? step_size_did : min(step_size_temp, max_step_size);
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
    dy = (*derivs)(x, y);
    SPDLOG_DEBUG_FILE("RKQC reaching: ");
    SPDLOG_DEBUG_FILE("  X = {:+9.8e};", x);
    SPDLOG_DEBUG_FILE("  dx = {:+9.8e};", step_size_did);
    for (int i = 0; i < DOF; i++) {
        SPDLOG_DEBUG_FILE(" Y[{0}] = {1:+9.8e}, dYdX = {2:+9.8e};", i, y[i], dy[i]);
    }
}

RungeKutta::STATUS RungeKutta::Solve(REAL step_start, REAL eps) {
    // * Initialize all relevant quantities
    SPDLOG_INFO_FILE("Start to solving the ODE with RungeKutta Method.");
    INIT();

    // * The initial point
    double x = _X[0];
    VD y = _Y[0];
    VD dydx = _dYdX[0];
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
        RKQC_SingleStep(x, y, dydx, step_size, eps, Y_Scale, step_size_did, step_size_next);
        _X.push_back(x);
        _Y.push_back(y);
        _Yeq.push_back(derivs->Yeq(x));
        _dYdX.push_back(dydx);

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
