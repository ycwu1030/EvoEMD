#include "EvoEMD/RungeKutta.h"

#include <iostream>

#include "EvoEMD/EMD.h"

using namespace std;
using namespace EvoEMD;
class ODE_TEST : public ODE_FUNCS {
public:
    ODE_TEST() : ODE_FUNCS() { Set_DOF(2); };
    ~ODE_TEST(){};

    VD dYdX(double x, VD y) { return {y[0] * (2 - y[1]), y[1] * (y[0] - 1)}; }
};

class FO_TEST : public ODE_FUNCS {
public:
    FO_TEST() : ODE_FUNCS() {
        Set_DOF(1);
        M = 100.0;
        g = 2.0;
        lambda = 1e5;
    }
    ~FO_TEST(){};

    VD dYdX(double x, VD y) {
        double yeq = Yield_Eq(M / x, M, g);
        return {-lambda / x / x * (y[0] * y[0] - yeq * yeq)};
    }

private:
    double M;
    double g;
    double lambda;
};

int main(int argc, char const *argv[]) {
    VD BOUND = {1, 2.7};
    ODE_TEST funcs;
    funcs.Set_X_BEGIN(0);
    funcs.Set_X_END(10);
    funcs.Set_BOUNDARY_CONDITION(BOUND);

    RungeKutta solver(&funcs);
    solver.Solve(0.01);
    solver.Dump_Solution("RK_test.dat");

    VD FOI = {Yield_Eq(100 / 1, 100, 2)};
    FO_TEST fo;
    fo.Set_X_BEGIN(1);
    fo.Set_X_END(100);
    fo.Set_BOUNDARY_CONDITION(FOI);

    RungeKutta fosolver(&fo);
    fosolver.Solve(0.01);
    fosolver.Dump_Solution("FO_test.dat");
    return 0;
}
