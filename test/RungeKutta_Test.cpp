#include "RungeKutta.h"
#include <iostream>

using namespace std;

class ODE_TEST: public ODE_FUNCS
{
public:
    ODE_TEST():ODE_FUNCS(){Set_DOF(2);};
    ~ODE_TEST(){};

    VD dYdX(double x, VD y){
        return {y[0]*(2-y[1]),y[1]*(y[0]-1)};
    }
};

int main(int argc, char const *argv[])
{
    VD BOUND = {1,2.7};
    ODE_TEST funcs;
    funcs.Set_X_BEGIN(0);
    funcs.Set_X_END(10);
    funcs.Set_BOUNDARY_CONDITION(BOUND);

    RungeKutta solver(&funcs);
    solver.Solve(0.01);
    solver.Dump_Solution("RK_test.dat");
    return 0;
}
