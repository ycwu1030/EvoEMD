#include <cmath>
#include <iostream>

#include "EvoEMD/EvoEMD.h"

using namespace std;
using namespace EvoEMD;
class HODE : public ODE_FUNCS {
public:
    HODE()
        : BR(1.0),
          ODE_FUNCS()  //, X_BEGIN(log(1.1e-16)), X_END(log(1e8)), BOUNDARY_CONDITION(3)
    {
        Set_DOF(2);
        X_BEGIN = log(1e-8);
        X_END = log(10);
        Y_BEGIN[0] = 1e-2;
        Y_BEGIN[1] = 1e-8;
    }
    ~HODE(){};

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
    REAL BR;
};

int main(int argc, char const *argv[]) {
    HODE hub;
    RungeKutta solver(&hub);
    solver.Solve(1e-3, 1e-4);
    solver.Dump_Solution("Hubble_rho_Coupled.dat");
    return 0;
}
