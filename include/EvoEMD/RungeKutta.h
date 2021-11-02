#ifndef __RUNGEKUTTA_H__
#define __RUNGEKUTTA_H__

#include <string>

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

/*
 * The ODE System
 */
class ODE_FUNCS {
protected:
    int DOF;
    REAL X_BEGIN;
    REAL X_END;
    VD Y_BEGIN;
    VD Delta_Y_Ratio_BEGIN;

public:
    ODE_FUNCS() : DOF(-1){};
    virtual ~ODE_FUNCS(){};

    void Set_DOF(int DOF) {
        this->DOF = DOF;
        Y_BEGIN = VD(DOF, 0);
        Delta_Y_Ratio_BEGIN = VD(DOF, 0);
    }
    void Set_BOUNDARIES(REAL X_BEGIN, REAL X_END, VD boundary, bool relative_to_equilibrium = false);

    int Get_DOF() const { return DOF; }
    REAL Get_X_BEGIN() const { return X_BEGIN; }
    REAL Get_X_END() const { return X_END; }
    VD Get_Y_BEGIN() const { return Y_BEGIN; }
    VD Get_Delta_Y_Ratio_BEGIN() const { return Delta_Y_Ratio_BEGIN; }
    VD operator()(REAL x, const VD &y, const VD &delta_y_ratio) { return dYdX(x, y, delta_y_ratio); }
    virtual VD dYdX(REAL x, const VD &y, const VD &delta_y_ratio) = 0;  // delta_y_ratio = 1 - Y/Yeq;
    virtual VD Yeq(REAL x) = 0;
    virtual VB Is_Thermalized() = 0;  // Check whether the particle is starting in thermalization or not
    virtual VB Can_be_Negative() = 0;
    virtual VB Should_be_Thermalized(REAL x, const VD &y, const VD &delta_y_ratio) = 0;
};

struct RK_Point {
    int DOF;
    REAL X;
    VD Y;
    VD dYdX;
    VD Yeq;             // Equilibrium value of Y
    VD Delta_Y_Ratio;   // 1 - Y/Yeq;
    VB Thermal_Status;  // Whether Y is following Yeq;
    RK_Point(int dof)
        : DOF(dof), Y(dof, 0), dYdX(dof, 0), Yeq(dof, 0), Delta_Y_Ratio(dof, 0), Thermal_Status(dof, false) {}
};
/*
 * RungeKutta method solving ODE:
 *  dY_i/dx = ODE_i(x,Y_j,...)
 */
class RungeKutta {
public:
    typedef enum { SUCCESS = 1, TOOMANYSTEP = -1 } STATUS;
    RungeKutta();
    RungeKutta(ODE_FUNCS *ode);
    ~RungeKutta() = default;

    /*
     * @brief Setting the ODE system
     * @param ode a pointer to the ODE object.
     */
    void Set_ODE(ODE_FUNCS *ode);

    /*
     * @brief Solve the ODE system
     * @param step_start, the initial step size
     * @param eps, tolarance in solving the ODE function
     */
    STATUS Solve(REAL step_start, REAL eps = 1e-6);

    /*
     * @brief dump the solution into file
     * @param filename, the file name
     */
    void Dump_Solution(std::string filename);

    VD Get_Solution_X() const { return _X; }
    VVD Get_Solution_Y() const { return _Y; }
    VD Get_Solution_Y_END() const { return _Y.back(); }

private:
    // * Following are parameters controling the adaptive step size
    REAL SAFETY;
    REAL POW_GROW;
    REAL POW_SHRINK;
    REAL TINY_REL_THRESHOLD;
    int MAXSTEPS;

    bool SOLVED;             // * Store the status whether the ODE functions is solved
    int DOF;                 // * The degree of freedon of the ODE functions
    REAL X_BEGIN;            // * The starting point in x (the argument)
    REAL X_END;              // * The ending point in x (the argument)
    VD Y_BEGIN;              // * The starting point of y
    VD Delta_Y_Ratio_BEGIN;  // * 1-y/yeq at the begin;

    VD _X;               // * The vector storing the points in x
    VVD _Y;              // * The vector storing the points in y, len(_Y) = len(_X) and the second dimension is DOF
    VVD _Delta_Y_Ratio;  // * The vector storing the value 1 - y/yeq
    VVD _Yeq;            // * The vector storing the equilibrium value of Y;
    VVD _dYdX;           // * Similar to _Y, but storing dy/dx

    VD ZERO_THRESHOLD;  // * The threshold to treat the corresponding value as 0;

    ODE_FUNCS *
        derivs;  // * The ODE function, dy/dx = derivs(x,y,xxxx), we don't own it, we just have one pointer pointing it.

    /*
     * @brief Clear the system.
     */
    void Clear();

    /*
     * @brief Initialize the system just before we solve the ODE
     *          We can change anything we like, before calling Solve()
     */
    void INIT();

    /**
     * @brief  4th-order RungeKutta Method
     * @note
     * @param  &p_cur: starting point
     * @param  &p_next: ending point
     * @param  step_size: used step size
     * @param  &should_follow_equilibrium: whether any component should be treated as following equilibrium
     * @param  &follow_equilibrium: whether the component is following equilibrium in this step
     * @retval whether any component that is not true in should_follow_equilibrium becomes equilibrium
     */
    bool RK4_SingleStep(const RK_Point &p_cur, RK_Point &p_next, const REAL step_size, VB &should_follow_equilibrium,
                        VB &follow_equilibrium);

    /**
     * @brief The Runge-Kutta method taking one step forward
     *        Based on 4th order Runge-Kutta and adding quality control (adaptive stepsize),
     *        such that we can kind of achieve 5th order accuracy
     * @param x, input: current value of x, output: next value of x
     * @param y, input: current value of y, output: next value of y
     * @param yeq, input: current value of yeq, output: next value of yeq
     * @param dydx, input: current value of dy/dx, output: next value of dy/dx
     * @param delta_y_ratio, input: current value of 1-y/yeq, output: next value of 1-y/yeq
     * @param thermal_status, input: current status of thermalization, output: next status of thermalization
     * @param step_size_guess, input: initial guess of current step size
     * @param eps, input: tolerance
     * @param Y_Scale, input: the scale in y to set the error (it is possible in different direction of y, the scale
     * is quite different)
     * @param step_size_did, output: actual step size we take
     * @param step_size_further, output: based on what we did, the prospect step size for next step
     */
    bool RKQC_SingleStep(REAL &x, VD &y, VD &yeq, VD &dydx, VD &delta_y_ratio, VB &thermal_status,
                         const REAL step_size_guess, const REAL eps, const VD &Y_Scale, REAL &step_size_did,
                         REAL &step_size_further);
    bool RKQC_SingleStep(const RK_Point &p_cur, const REAL step_size_guess, const REAL eps, RK_Point &p_next,
                         REAL &step_size_did, REAL &step_size_further);
};
}  // namespace EvoEMD
#endif  //__RUNGEKUTTA_H__
