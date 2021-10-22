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
    VD BOUNDARY_CONDITION;

public:
    ODE_FUNCS(){};
    virtual ~ODE_FUNCS(){};

    void Set_DOF(int dof) { DOF = dof; }
    void Set_X_BEGIN(REAL xi) { X_BEGIN = xi; }
    void Set_X_END(REAL xf) { X_END = xf; }
    void Set_BOUNDARY_CONDITION(VD boundary) { BOUNDARY_CONDITION = boundary; }

    int Get_DOF() const { return DOF; }
    REAL Get_X_BEGIN() const { return X_BEGIN; }
    REAL Get_X_END() const { return X_END; }
    VD Get_BOUNDARY_CONDITION() const { return BOUNDARY_CONDITION; }
    VD operator()(REAL x, VD y) { return dYdX(x, y); }
    virtual VD dYdX(REAL x, VD y) = 0;
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

    bool SOLVED;   // * Store the status whether the ODE functions is solved
    int DOF;       // * The degree of freedon of the ODE functions
    REAL X_BEGIN;  // * The starting point in x (the argument)
    REAL X_END;    // * The ending point in x (the argument)

    VD _X;      // * The vector storing the points in x
    VVD _Y;     // * The vector storing the points in y, len(_Y) = len(_X) and the second dimension is DOF
    VVD _dYdX;  // * Similar to _Y, but storing dy/dx

    VD BOUNDARY_AT_BEGIN;  // * The boundary condition at X_BEGIN
    VD ZERO_THRESHOLD;     // * The threshold to treat the corresponding value as 0;

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
     * @brief The 4th order Runge-Kutta method taking one step forward
     * @param x_cur, current x value
     * @param y_cur, current y value
     * @param dy_cur, current dy/dx value
     * @param step_size, step size, we would like to go forward
     * @param dy_next, the dy at next step (output)
     */
    void RK4_SingleStep(const REAL x_cur, const VD &y_cur, const VD &dy_cur, const REAL step_size, VD &y_next);

    /**
     * @brief The Runge-Kutta method taking one step forward
     *        Based on 4th order Runge-Kutta and adding quality control (adaptive stepsize),
     *        such that we can kind of achieve 5th order accuracy
     * @param x, input: current value of x, output: next value of x
     * @param y, input: current value of y, output: next value of y
     * @param dy, input: current value of dy, output: next value of dy
     * @param step_size_guess, input: initial guess of current step size
     * @param eps, input: tolerance
     * @param Y_Scale, input: the scale in y to set the error (it is possible in different direction of y, the scale is
     * quite different)
     * @param step_size_did, output: actual step size we take
     * @param step_size_further, output: based on what we did, the prospect step size for next step
     */
    void RKQC_SingleStep(REAL &x, VD &y, VD &dy, const REAL step_size_guess, const REAL eps, const VD &Y_Scale,
                         REAL &step_size_did, REAL &step_size_further);
};
}  // namespace EvoEMD
#endif  //__RUNGEKUTTA_H__
