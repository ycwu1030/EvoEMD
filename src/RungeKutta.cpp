#include "RungeKutta.h"
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

RungeKutta::RungeKutta():DOF(0),derivs(nullptr),SOLVED(false),SAFETY(0.9),POW_GROW(-0.2),POW_SHRINK(-0.25),TINY(1e-30),MAXSTEPS(10000)
{ }

RungeKutta::RungeKutta(ODE_FUNCS *ode):SOLVED(false),SAFETY(0.9),POW_GROW(-0.2),POW_SHRINK(-0.25),TINY(1e-30),MAXSTEPS(10000)
{
    Set_ODE(ode);
}

void RungeKutta::Set_ODE(ODE_FUNCS *ode)
{
    derivs = ode;
    SOLVED = false;
}

void RungeKutta::Clear()
{
    _X.clear();
    _Y.clear();
    _dYdX.clear();
}

void RungeKutta::INIT()
{
    // * Getting information from the ODE system
    DOF = derivs->Get_DOF();
    X_BEGIN = derivs->Get_X_BEGIN();
    X_END = derivs->Get_X_END();
    BOUNDARY_AT_BEGIN = derivs->Get_BOUNDARY_CONDITION();
    
    // * Clear whatever we have in our system
    Clear();

    // * Initialize the system at X_BEGIN:
    _X.push_back(X_BEGIN);
    _Y.push_back(BOUNDARY_AT_BEGIN);
    _dYdX.push_back((*derivs)(X_BEGIN,BOUNDARY_AT_BEGIN));
}

void RungeKutta::RK4_SingleStep(const REAL x_cur, const VD &y_cur, const VD &dy_cur, const REAL step_size, VD &y_next)
{
    REAL half_step = step_size/2.0;

    // * 1. Using dx*dy/dx
    VD dY_Step1 = step_size*dy_cur;

    // * 2. Using dx/2 and dY_Step1/2
    VD dY_Step2 = step_size*(*derivs)(x_cur + half_step, y_cur + dY_Step1/2.0);

    // * 3. Using dx/2 and dY_Step2/2
    VD dY_Step3 = step_size*(*derivs)(x_cur + half_step, y_cur + dY_Step2/2.0);

    // * 4. Using dx and dY_Step3
    VD dY_Step4 = step_size*(*derivs)(x_cur + step_size, y_cur + dY_Step3);

    // * Combine above 4 steps:
    y_next = y_cur + dY_Step1/6.0 + dY_Step2/3.0 + dY_Step3/3.0 + dY_Step4/6.0;
}

void RungeKutta::RKQC_SingleStep(REAL &x, VD &y, VD &dy, const REAL step_size_guess, const REAL eps, const VD &Y_Scale, REAL &step_size_did, REAL &step_size_further)
{
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
    REAL min_step_size = 1e-5*step_size_guess;
    REAL max_step_size;

    while (true)
    {
        // * Take two half steps
        half_step_size = step_size/2.0;
        RK4_SingleStep(x_cache,y_cache,dy_cache,half_step_size,y_temp);
        x = x_cache + half_step_size;
        dy = (*derivs)(x,y_temp);
        RK4_SingleStep(x,y_temp,dy,half_step_size,y);

        // * Take one full step
        RK4_SingleStep(x_cache,y_cache,dy_cache,step_size,y_temp);
        x = x_cache + step_size;

        // * Check the difference between above two methods
        Delta_y = y - y_temp;
        error_temp = fabs(Delta_y/Y_Scale);
        error_max = *max_element(error_temp.begin(),error_temp.end());
        error_max /= eps;
        if (error_max <= 1.0)
        {
            // * The error is acceptable, we will proceed, and even try to enlarge the step size
            step_size_did = step_size;
            step_size_temp = SAFETY*step_size*exp(POW_GROW*log(error_max));
            max_step_size = 4*step_size_did;
            step_size_further = min(step_size_temp,max_step_size);
            break;
        }
        // * Otherwise, the error is larger than our tolarance, we need to shrink the step size to improve the error
        step_size_temp = SAFETY*step_size*exp(POW_SHRINK*log(error_max));
        if (fabs(step_size_temp) < fabs(min_step_size))
        {
            // * If after shrink, the step size is too small, then we don't shrink it any further, and accept current possible large error.
            step_size_did = step_size;
            step_size_further = step_size;
            break;
        }
        step_size = step_size_temp; // * Change the step size, do the RK4 again.
    }
    y = y + Delta_y/15;
    dy = (*derivs)(x,y);
}

RungeKutta::STATUS RungeKutta::Solve(REAL step_start, REAL eps)
{
    // * Initialize all relevant quantities
    INIT();

    // * The initial point
    double x = _X[0];
    VD y = _Y[0];
    VD dydx = _dYdX[0];
    VD Y_Scale(DOF);

    double step_size = (X_END > X_BEGIN)?fabs(step_start):-fabs(step_start);
    double step_size_did;
    double step_size_next;

    for (int nstp = 0; nstp < MAXSTEPS; nstp++)
    {
        Y_Scale = fabs(y) + fabs(dydx*step_size);
        for (int i = 0; i < DOF; i++)
        {
            Y_Scale[i] = min(1.0,Y_Scale[i]);
        }
        
        // * Check whether the step size is too large that we already pass the end point
        if ( (step_size > 0 && x+step_size > X_END) || (step_size < 0 && x + step_size < X_END) )
        {
            step_size = X_END - x;
        }

        // * One step forward
        RKQC_SingleStep(x,y,dydx,step_size,eps,Y_Scale,step_size_did,step_size_next);
        _X.push_back(x);
        _Y.push_back(y);
        _dYdX.push_back(dydx);
        
        // * Check if we reach the end point
        if ( (x-X_END)*(X_END-X_BEGIN) >= 0 )
        {
            SOLVED = true;
            return SUCCESS;
        }
        
        // * If not reaching the end, we adapt the stepsize, but need to control it either not too large or too small.
        step_size = step_size_next;
        if (fabs(step_size_next) > 1e-1*fabs(X_END-X_BEGIN))
        {
            step_size = 1e-1*(X_END-X_BEGIN);
        }
        if (fabs(step_size_next) < 1e-5*fabs(X_END-X_BEGIN))
        {
            step_size = 1e-5*(X_END-X_BEGIN);
        }
    }
    SOLVED = true;
    return TOOMANYSTEP;
}

void RungeKutta::Dump_Solution(string filename)
{
    ofstream output(filename.c_str());
    output<<"x\t";
    for (size_t i = 0; i < DOF; i++)
    {
        output<<"y_"<<i<<"\t";
    }
    for (size_t i = 0; i < DOF; i++)
    {
        output<<"dy_"<<i<<"/dx"<<"\t";
    }
    output<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        output<<_X[i]<<"\t";
        for (size_t j = 0; j < DOF; j++)
        {
            output<<_Y[i][j]<<"\t";
        }
        for (size_t j = 0; j < DOF; j++)
        {
            output<<_dYdX[i][j]<<"\t";
        }
        output<<endl;
    }
    output.close();
}
