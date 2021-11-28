#ifndef _HUBBLE_EVO_H_
#define _HUBBLE_EVO_H_

#include <cmath>

#include "EvoEMD/Hubble.h"
#include "EvoEMD/ParameterBase.h"
#include "EvoEMD/RungeKutta.h"
#include "gsl/gsl_spline.h"

namespace EvoEMD {

class Hubble_Evolution : public ODE_FUNCS {
    // * Used to solve ODE of
    // * dY1/du = - u Y1/sqrt(u Y1 + Y2);
    // * dY2/du = f u^2 Y1/sqrt(u Y1 + Y2);
    // * Above functions are invariant under following scaling transformation
    // * Y1 -> k^3 Y1;
    // * Y2 -> k^4 Y2;
    // * u  -> k u;
    // * Hence, the initial conditions can be implemented as Y1/u^3 = xxx and Y2/u^4 = xxx.
    // * Then one can start solving the equation with any value of u;
    // * In practice we will solve above equations in terms of x = ln(u)
    // * dY1/dx = - exp(2x) Y1/sqrt(exp(x)Y1 + Y2)
    // * dY2/dx = f exp(3x) Y1/sqrt(exp(x)Y1 + Y2)
public:
    Hubble_Evolution();
    ~Hubble_Evolution(){};

    void Solve(REAL Uini, REAL Y1ini, REAL Y2ini, REAL BR = 1);

    VD Get_Solution_Y1() const { return Y1; }
    VD Get_Solution_Y2() const { return Y2; }
    VD Get_Solution_U() const { return U; }

    virtual VD dYdX(REAL x, const VD &y, const VD &delta_y_ratio) {
        VD res(DOF);
        res[0] = -exp(2 * x) * y[0] / sqrt(exp(x) * y[0] + y[1]);
        res[1] = BR * exp(3 * x) * y[0] / sqrt(exp(x) * y[0] + y[1]);
        return res;
    }
    virtual VD Yeq(REAL x) { return VD(DOF, 0); }
    virtual VB Is_Thermalized() { return VB(DOF, false); }
    virtual VB Can_be_Negative() { return VB(DOF, false); }
    virtual VB Should_be_Thermalized(REAL x, const VD &y, const VD &delta_y_ratio) { return VB(DOF, false); }

private:
    VD X;  // = log(U)
    VD U;
    VD Y1;
    VD Y2;
    REAL BR;

    RungeKutta rk;
};

class Hubble_BE : public Parameter_Base, public Hubble_For_Single_Period {
private:
    VD List_T;
    VD List_U;  // U = k*a
    VD List_H;
    VD List_rhoR;
    VD List_rhoM;

    gsl_interp_accel *acc_Us;
    gsl_interp_accel *acc_Hs;
    gsl_spline *spline_Us;
    gsl_spline *spline_Hs;

    Hubble_BE();
    ~Hubble_BE(){};

    Hubble_Evolution HE;
    REAL Calc_Tr(VD U, VD Y1, VD Y2, REAL krhoR);
    void Clean();

public:
    static Hubble_BE &Get_Hubble_Calculator();
    virtual REAL Get_Hubble_at_T(const REAL T) override;
    virtual REAL Get_dlna_dlnT_at_T(const REAL T) override;
    virtual void Update_Value(REAL input) override;
    VD Get_T_List();
    VD Get_U_List();
    VD Get_H_List();
    VD Get_rhoR_List();
    VD Get_rhoM_List();
};

}  // namespace EvoEMD

#endif  //_HUBBLE_EVO_H_
