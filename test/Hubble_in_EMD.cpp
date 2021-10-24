#include <cmath>
#include <iostream>

#include "EvoEMD/EMD.h"
#include "EvoEMD/RungeKutta.h"

using namespace std;
using namespace EvoEMD;
class HODE : public ODE_FUNCS {
public:
    HODE()
        : ODE_FUNCS()  //, X_BEGIN(log(1.1e-16)), X_END(log(1e8)), BOUNDARY_CONDITION(3)
    {
        Set_DOF(3);
        c1 = sqrt(3) * 5.0 / 2.0 * M_PI * sqrt(106.5) / 3.0 / sqrt(10);
        cs = c1 * pow(2 * M_PI * M_PI / 45 * 106.5, 1.0 / 3.0);
        X_BEGIN = log(1e-30);
        X_END = log(1e20);
        BOUNDARY_CONDITION = VD(3);
        BOUNDARY_CONDITION[1] = 1790;
        BOUNDARY_CONDITION[2] = 897;
        BOUNDARY_CONDITION[0] = BOUNDARY_CONDITION[1] / 1.1e-16;
    }
    ~HODE(){};

    virtual VD dYdX(REAL x, VD y) {
        VD res(DOF);
        REAL A = exp(x);
        res[0] = -c1 * A * A * y[0] / sqrt(A * y[0] + y[1]);
        res[1] = c1 * A * A * A * y[0] / sqrt(A * y[0] + y[1]);
        res[2] = cs * A * A * A * y[0] / pow(y[2], 1.0 / 3.0) / sqrt(A * y[0] + y[1]);
        return res;
    }

private:
    REAL c1;
    REAL cs;
};

int main(int argc, char const *argv[]) {
    HODE hub;
    RungeKutta solver(&hub);
    solver.Solve(0.01);
    solver.Dump_Solution("Hubble_Coupled.dat");
    return 0;
}
