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

#define LOGGING_RK_POINT(rkp)                                                                                         \
    SPDLOG_DEBUG_FILE("  X = {:+9.8e}", rkp.X);                                                                       \
    for (int irkp = 0; irkp < rkp.DOF; ++irkp) {                                                                      \
        SPDLOG_DEBUG_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e}", irkp, rkp.Y[irkp], rkp.dYdX[irkp]);          \
        SPDLOG_DEBUG_FILE("      Yeq[{0}] = {1:+9.8e}, dR[{0}] = {2:+9.8e}, Thermal[{0}] = {3}", irkp, rkp.Yeq[irkp], \
                          rkp.Delta_Y_Ratio[irkp], rkp.Thermal_Status[irkp]);                                         \
    }

struct RK_INTER {
    VD dY_Step1;
    VD dY_Step2;
    VD dY_Step3;
    VD dY_Step4;

    VD Y1;
    VD dYR1;
    VD Y2;
    VD dYR2;
    VD Y3;
    VD dYR3;

    VD Yeq_Mid;
    VD Yeq_End;
    RK_INTER(int dof)
        : dY_Step1(dof, 0),
          dY_Step2(dof, 0),
          dY_Step3(dof, 0),
          dY_Step4(dof, 0),
          Y1(dof, 0),
          dYR1(dof, 0),
          Y2(dof, 0),
          dYR2(dof, 0),
          Y3(dof, 0),
          dYR3(dof, 0),
          Yeq_Mid(dof, 0),
          Yeq_End(dof, 0) {}
};

#define LOGGING_RK_INTER(rki)                                                                           \
    SPDLOG_DEBUG_FILE("The intermediate status: ");                                                     \
    for (int irki = 0; irki < rki.dY_Step1.size(); irki++) {                                            \
        SPDLOG_DEBUG_FILE("For Component-{}", irki);                                                    \
        SPDLOG_DEBUG_FILE("dY1={:+9.8e}, dY2={:+9.8e}, dY3={:+9.8e}, dY4={:+9.8e}", rki.dY_Step1[irki], \
                          rki.dY_Step2[irki], rki.dY_Step3[irki], rki.dY_Step4[irki]);                  \
    }

/**
 * @brief  Checking whether the component is thermalized during this step.
 * @note
 * @param  &p_cur: starting point of this step
 * @param  &p_inter: intermediate states of this step
 * @param  &should_follow_equilibrium: whether the component is already keep thermalized within this step
 * @retval The vector<bool> indicates whether this component is thermalized but not considered in follow_equilibrium
 */
VB WHETHER_THERMALIZED(const RK_Point &p_cur, const RK_INTER &p_inter, const VB &should_follow_equilibrium) {
    VB res(p_cur.DOF, false);
    for (int i = 0; i < p_cur.DOF; i++) {
        if (!should_follow_equilibrium[i]) {
            // * Only check when we didn't considered it as following equilibrium during this step
            // * Then we have two cases:
            // * - it starts with thermalized: freeze-out like
            // * - it starts with un-thermalized: freeze-in like
            if (p_cur.Thermal_Status[i]) {
                // * Freeze-out like case;
                // * In most cases (maybe all), dY_Step1 will be zero
                // * So we start checking from dY_Step2
                if (signbit(p_inter.dY_Step2[i] / p_inter.dY_Step3[i]) &&
                    signbit(p_inter.dY_Step3[i] / p_inter.dY_Step4[i])) {
                    REAL Y_temp = p_cur.Y[i] + p_inter.dY_Step1[i] / 6.0 + p_inter.dY_Step2[i] / 3.0 +
                                  p_inter.dY_Step3[i] / 3.0 + p_inter.dY_Step4[i] / 6.0;
                    if (Y_temp < p_inter.Yeq_End[i]) {
                        res[i] = true;
                    }
                }
            } else {
                // * Freeze-in like case;
                // * We have two cases:
                // * - starting below equilibrium
                // * - starting above equilibrium
                if (p_cur.Delta_Y_Ratio[i] > 0) {
                    // * Below equilibrium
                    if (p_inter.Y1[i] > p_inter.Yeq_Mid[i] &&
                        p_inter.Y1[i] + p_inter.dY_Step2[i] / 2.0 < p_inter.Yeq_End[i]) {
                        res[i] = true;
                    }
                } else {
                    // * Above equilibrium
                    if (p_inter.Y1[i] < p_inter.Yeq_Mid[i] &&
                        p_inter.Y1[i] + p_inter.dY_Step2[i] / 2.0 > p_inter.Yeq_End[i]) {
                        res[i] = true;
                    }
                }
            }
        }
    }
    return res;
}

bool RungeKutta::RK4_SingleStep(const RK_Point &p_cur, RK_Point &p_next, const REAL step_size,
                                VB &should_follow_equilibrium, VB &follow_equilibrium) {
    SPDLOG_INFO_FILE("4th Order RK with stepsize = {:+9.8e}.", step_size);
    SPDLOG_INFO_FILE("Starting at: ");
    LOGGING_RK_POINT(p_cur);
    RK_INTER inter(DOF);
    bool comp;
    VB should_be_thermalized = derivs->Should_be_Thermalized(p_cur.X, p_cur.Y, p_cur.Delta_Y_Ratio);

    // * 4th order RK proceeds as:
    // * x0 -> x1 -> x2
    // * x1 = x0 + dx/2
    // * x2 = x0 + dx
    // * Yeq_Mid = Yeq(x1)
    // * Yeq_End = Yeq(x2)
    // * 1st: dy1 = dx * dydx(x0, y0)
    // *      y1 = y0 + dy1/2;   y1 can be checked against Yeq_Mid;
    // * 2nd: dy2 = dx * dydx(x1, y1)
    // *      y2 = y0 + dy2/2;   y2p = y1 + dy2/2 can be checked against Yeq_End;
    // * 3rd: dy3 = dx * dydx(x1, y2)
    // *      y3 = y0 + dy3;
    // * 4th: dy4 = dx * dydx(x2, y3)
    // * x0 -> x2, y0 -> y0 + dy1/6 + dy2/3 + dy3/3 + dy4/6;
    REAL dx = step_size;
    REAL dx_half = step_size / 2.0;
    REAL x1 = p_cur.X + dx_half;
    REAL x2 = p_cur.X + dx;
    inter.Yeq_Mid = derivs->Yeq(x1);
    inter.Yeq_End = derivs->Yeq(x2);

    VB can_be_negative = derivs->Can_be_Negative();
    REAL Y_temp;
    while (true) {
        for (int i = 0; i < DOF; i++) {
            SPDLOG_DEBUG_FILE("Component-{}, should be thermalized? {}", i, should_be_thermalized[i]);
        }
        // * 1. Using dx*dy/dx
        inter.dY_Step1 = dx * p_cur.dYdX;
        for (int i = 0; i < DOF; i++) {
            if (should_be_thermalized[i]) {
                // * This component will follow the thermal equilibrium during current step
                inter.Y1[i] = inter.Yeq_Mid[i];
                inter.dYR1[i] = 0;
            } else {
                Y_temp = p_cur.Y[i] + inter.dY_Step1[i] / 2.0;
                inter.Y1[i] = (Y_temp < 0 && (!can_be_negative[i])) ? 0 : Y_temp;
                inter.dYR1[i] = 1.0 - inter.Y1[i] / inter.Yeq_Mid[i];
            }
        }

        // * 2. Using dx/2 and dY_Step1/2
        inter.dY_Step2 = dx * (*derivs)(x1, inter.Y1, inter.dYR1);
        for (int i = 0; i < DOF; i++) {
            if (should_be_thermalized[i]) {
                // * This component will follow the thermal equilibrium during current step
                inter.Y2[i] = inter.Yeq_Mid[i];
                inter.dYR2[i] = 0;
            } else {
                Y_temp = p_cur.Y[i] + inter.dY_Step2[i] / 2.0;
                inter.Y2[i] = (Y_temp < 0 && (!can_be_negative[i])) ? 0 : Y_temp;
                inter.dYR2[i] = 1.0 - inter.Y2[i] / inter.Yeq_Mid[i];
            }
        }

        // * 3. Using dx/2 and dY_Step2/2
        inter.dY_Step3 = dx * (*derivs)(x1, inter.Y2, inter.dYR2);
        for (int i = 0; i < DOF; i++) {
            if (should_be_thermalized[i]) {
                inter.Y3[i] = inter.Yeq_End[i];
                inter.dYR3[i] = 0;
            } else {
                Y_temp = p_cur.Y[i] + inter.dY_Step3[i];
                inter.Y3[i] = (Y_temp < 0 && (!can_be_negative[i])) ? 0 : Y_temp;
                inter.dYR3[i] = 1.0 - inter.Y3[i] / inter.Yeq_End[i];
            }
        }

        // * 4. Using dx and dY_Step3
        inter.dY_Step4 = dx * (*derivs)(x2, inter.Y3, inter.dYR3);

        LOGGING_RK_INTER(inter);
        // * Combine above 4 steps:
        bool match_thermal = true;
        p_next.DOF = DOF;
        p_next.X = x2;
        p_next.Yeq = inter.Yeq_End;
        for (int i = 0; i < DOF; i++) {
            if (should_be_thermalized[i]) {
                p_next.Y[i] = p_next.Yeq[i];
                p_next.Delta_Y_Ratio[i] = 0;
                p_next.Thermal_Status[i] = true;
                comp = true;
            } else {
                Y_temp = p_cur.Y[i] + inter.dY_Step1[i] / 6.0 + inter.dY_Step2[i] / 3.0 + inter.dY_Step3[i] / 3.0 +
                         inter.dY_Step4[i] / 6.0;
                if (p_cur.Thermal_Status[i]) {
                    if (Y_temp < p_next.Yeq[i]) {
                        match_thermal = false;
                        should_be_thermalized[i] = true;
                    }
                }
                p_next.Y[i] = (Y_temp < 0 && (!can_be_negative[i])) ? 0 : Y_temp;
                p_next.Delta_Y_Ratio[i] = 1.0 - p_next.Y[i] / p_next.Yeq[i];
                p_next.Thermal_Status[i] = false;
            }
        }
        if (match_thermal) break;
    }

    p_next.dYdX = derivs->dYdX(p_next.X, p_next.Y, p_next.Delta_Y_Ratio);
    follow_equilibrium = should_be_thermalized;
    return comp;
}

bool RungeKutta::RKQC_SingleStep(const RK_Point &p_cur, const REAL step_size_guess, const REAL eps, RK_Point &p_next,
                                 REAL &step_size_did, REAL &step_size_further) {
    SPDLOG_INFO_FILE("RKQC: ");
    SPDLOG_INFO_FILE("Starting with: ");
    SPDLOG_INFO_FILE("Guess step size: dx = {:+9.8e},", step_size_guess);
    LOGGING_RK_POINT(p_cur);

    // *                  p_temp
    // *              a /        \ b
    // *       (p_cache)           p_end
    // * p_cur (p_cache) --------  p_next
    // *                    c

    // * Cache the initial points
    RK_Point p_cache = p_cur;
    RK_Point p_temp(DOF);
    RK_Point p_end(DOF);
    REAL step_size = step_size_guess;
    REAL step_size_temp;
    REAL half_step_size;

    VD Delta_y(DOF, 0);
    VD error_temp(DOF, 0);

    REAL error_max = 0;
    REAL min_step_size = 1e-5 * step_size_guess;
    SPDLOG_DEBUG_FILE("Minimal allowed step size is {:+9.8e};", min_step_size);
    REAL max_step_size;
    int trials = 0;
    VB should_follow_equilibrium(DOF, false);
    VB follow_equilibrium_two_step_first(DOF, false);
    VB follow_equilibrium_two_step_second(DOF, false);
    VB follow_equilibrium_single_step(DOF, false);
    bool follow_equilibrium_any;
    bool new_component_thermalized;
    while (true) {
        ++trials;
        SPDLOG_DEBUG_FILE("{}-th Trial for RKQC with stepsize = {:+9.8e}", trials, step_size);
        // * Take two half steps
        half_step_size = step_size / 2.0;
        new_component_thermalized = RK4_SingleStep(p_cache, p_temp, half_step_size, should_follow_equilibrium,
                                                   follow_equilibrium_two_step_first);
        new_component_thermalized = RK4_SingleStep(p_temp, p_end, half_step_size, should_follow_equilibrium,
                                                   follow_equilibrium_two_step_second);

        SPDLOG_DEBUG_FILE("Take two half steps, reaching:");
        LOGGING_RK_POINT(p_end);
        for (int i = 0; i < DOF; i++) {
            if (p_end.Thermal_Status[i] && p_end.Y[i] < ZERO_THRESHOLD[i]) {
                SPDLOG_WARN_FILE("  REACH ZERO THRESHOLD AS Y[{0}] = {1:+9.8e} < {2:+9.8e};", i, p_end.Y[i],
                                 ZERO_THRESHOLD[i]);
                p_end.Y[i] = 0;
            }
        }

        // * Take one full step
        new_component_thermalized =
            RK4_SingleStep(p_cache, p_next, step_size, should_follow_equilibrium, follow_equilibrium_single_step);

        SPDLOG_DEBUG_FILE("Take one full step, reaching:");
        LOGGING_RK_POINT(p_next);
        for (int i = 0; i < DOF; i++) {
            if (p_next.Thermal_Status[i] && p_next.Y[i] < ZERO_THRESHOLD[i]) {
                SPDLOG_WARN_FILE("  REACH ZERO THRESHOLD AS Y[{0}] = {1:+9.8e} < {2:+9.8e};", i, p_next.Y[i],
                                 ZERO_THRESHOLD[i]);
                p_next.Y[i] = 0;
            }
        }

        // * Checking if there is any component is following equilibrium, and if that component drops too fast;
        // * Also checking if we are crossing the critical point
        follow_equilibrium_any = false;
        bool drop_fast_or_crossing_critical = false;
        for (int i = 0; i < DOF; i++) {
            if (follow_equilibrium_two_step_second[i] && follow_equilibrium_single_step[i]) {
                follow_equilibrium_any = true;
                // * This component follows equilibrium
                REAL drop_slop = p_next.Y[i] / p_cache.Y[i];
                if (drop_slop < 0.2) {
                    // * It drops too fast, we have to shrink the step to be safe.
                    drop_fast_or_crossing_critical = true;
                }
            } else if (!follow_equilibrium_two_step_second[i] && follow_equilibrium_single_step[i]) {
                // * In step b, this component is not following equilibrium
                drop_fast_or_crossing_critical = true;
            }
        }
        if (drop_fast_or_crossing_critical) {
            step_size_temp = step_size;
            step_size = step_size_temp / 2.0;
            SPDLOG_WARN_FILE(
                "Either when following the equilibrium, it drops too fast, or we probably cross the critical point, we "
                "need to shrink the step size to be safe: {:+9.8e} -> {:+9.8e}",
                step_size_temp, step_size);
            continue;
        }

        // * If we don't have drop too fast and crossing critical point problem, then we have to check the error

        // ** First calculate the error
        VD Y_Scale = fabs(p_cur.Y) + fabs(p_cur.dYdX * step_size);
        for (int i = 0; i < DOF; i++) {
            Y_Scale[i] = max(1e-50, Y_Scale[i]);
            SPDLOG_INFO_FILE("  Y[{0}] = {1:+9.8e}, dYdX[{0}] = {2:+9.8e}, YScale[{0}] = {3:+9.8e};", i, p_cur.Y[i],
                             p_cur.dYdX[i], Y_Scale[i]);
        }
        Delta_y = p_end.Y - p_next.Y;
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
            step_size_further = (follow_equilibrium_any) ? 1.5 * step_size_did : min(step_size_temp, max_step_size);
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
    p_next.Y = p_end.Y + Delta_y / 15.0;
    for (int i = 0; i < DOF; i++) {
        if (p_next.Thermal_Status[i]) {
            p_next.Delta_Y_Ratio[i] = 0;
        } else {
            p_next.Delta_Y_Ratio[i] = 1.0 - p_next.Y[i] / p_next.Yeq[i];
        }
    }
    p_next.dYdX = (*derivs)(p_next.X, p_next.Y, p_next.Delta_Y_Ratio);
    SPDLOG_DEBUG_FILE("RKQC reaching: ");
    LOGGING_RK_POINT(p_next);
    return follow_equilibrium_any;
}

RungeKutta::STATUS RungeKutta::Solve(REAL step_start, REAL eps) {
    // * Initialize all relevant quantities
    SPDLOG_INFO_FILE("Start to solving the ODE with RungeKutta Method.");
    INIT();

    // * The initial point
    RK_Point p_start(DOF);
    RK_Point p_end(DOF);
    p_start.X = _X[0];
    p_start.Y = _Y[0];
    p_start.Yeq = _Yeq[0];
    p_start.dYdX = _dYdX[0];
    p_start.Delta_Y_Ratio = _Delta_Y_Ratio[0];
    p_start.Thermal_Status = derivs->Is_Thermalized();

    double step_size = (X_END > X_BEGIN) ? std::fabs(step_start) : -std::fabs(step_start);
    double step_size_did;
    double step_size_next;

    for (int nstp = 0; nstp < MAXSTEPS; nstp++) {
        SPDLOG_INFO_FILE("Step-{}:", nstp + 1);
        LOGGING_RK_POINT(p_start);

        // * Check whether the step size is too large that we already pass the end point
        if ((step_size > 0 && p_start.X + step_size > X_END) || (step_size < 0 && p_start.X + step_size < X_END)) {
            step_size = X_END - p_start.X;
        }

        // * One step forward
        bool comp = RKQC_SingleStep(p_start, step_size, eps, p_end, step_size_did, step_size_next);

        _X.push_back(p_end.X);
        _Y.push_back(p_end.Y);
        _Yeq.push_back(p_end.Yeq);
        _dYdX.push_back(p_end.dYdX);
        _Delta_Y_Ratio.push_back(p_end.Delta_Y_Ratio);
        p_start = p_end;

        // * Check if we reach the end point
        if ((p_end.X - X_END) * (X_END - X_BEGIN) >= 0) {
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
