#include "LeptogenesisBE.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

LeptogenesisBE::LeptogenesisBE() : ODE_FUNCS(), EMDEvo(1e14, 1e-1)
{
    Set_DOF(R_Calculator.Get_DOF()); // N1, N2, lep-number, chi, S
    solver.Set_ODE(this);
}

void LeptogenesisBE::Set_Temperatures(REAL Ti, REAL Tr, REAL Tf, REAL TInflation)
{
    EMDEvo.Set_Temperature(Ti, Tr, Tf, TInflation);
}

void LeptogenesisBE::Set_Nu_Masses(REAL MNR1, REAL MNR2, REAL MNR3)
{
    R_Calculator.Set_Nu_Masses(MNR1, MNR2, MNR3);
}

void LeptogenesisBE::Set_Nu_MassOrdering(Nu_TypeI_SeeSaw::MassOrdering od)
{
    R_Calculator.Set_Nu_MassOrdering(od);
}

void LeptogenesisBE::Set_Dark_Sector_Masses(REAL MCHI, REAL MS)
{
    R_Calculator.Set_Dark_Sector_Masses(MCHI, MS);
}

void LeptogenesisBE::Set_LambdaX(REAL lamx)
{
    R_Calculator.Set_LambdaX(lamx);
}

VD LeptogenesisBE::dYdX(REAL x, VD y)
{
    switch (PeriodID)
    {
    case 0:
        return dYdX0(x, y);
    case 1:
        return dYdX1(x, y);
    case 2:
        return dYdX2(x, y);
    case 3:
        return dYdX3(x, y);
    default:
        return dYdX3(x, y);
    }
}

VD LeptogenesisBE::dYdX0(REAL x, VD y)
{
    REAL z = exp(x);
    VD res(DOF);
    REAL Temp = R_Calculator.Get_Temperature_from_Z(z);
    REAL Temp1 = R_Calculator.Get_Temperature_from_Z(1);
    res = z * z * R_Calculator.Ri(Temp, y) / Entropy_Density(Temp) / EMDEvo.Get_Hubble_at_T_Period(EMD::ERDE, Temp1);
    return res;
}

VD LeptogenesisBE::dYdX1(REAL x, VD y)
{
    REAL z = exp(x);
    VD res(DOF);
    REAL Temp = R_Calculator.Get_Temperature_from_Z(z);
    REAL Temp1 = R_Calculator.Get_Temperature_from_Z(1);
    res = z * sqrt(z) * R_Calculator.Ri(Temp, y) / Entropy_Density(Temp) / EMDEvo.Get_Hubble_at_T_Period(EMD::EMDE, Temp1);
    return res;
}

VD LeptogenesisBE::dYdX2(REAL x, VD y)
{
    REAL z = exp(x);
    VD res(DOF);
    REAL Temp = R_Calculator.Get_Temperature_from_Z(z);
    REAL Temp1 = R_Calculator.Get_Temperature_from_Z(1);
    REAL zr = R_Calculator.Get_Z_from_Temperature(EMDEvo.Get_Tr());
    REAL HEPatZ = EMDEvo.Get_Hubble_at_T_Period(EMD::EPE, Temp);
    REAL HRDati = EMDEvo.Get_Hubble_at_T_Period(EMD::ERDE, EMDEvo.Get_Ti());
    REAL HRDatr = EMDEvo.Get_Hubble_at_T_Period(EMD::RDE, EMDEvo.Get_Tr());
    REAL HEPat1 = EMDEvo.Get_Hubble_at_T_Period(EMD::EPE, Temp1);
    REAL satz = Entropy_Density(Temp);
    REAL sati = Entropy_Density(EMDEvo.Get_Ti());
    REAL entropy_term = 20.0 * EMDEvo.Get_Delta() / 3.0 * pow(HEPatZ / HRDati, 2) * HRDatr / HEPat1 * pow(z, 4) / zr * sati / satz * exp(-5.0 * HRDatr / HEPatZ / 3.0);
    res = 8.0 * z * z * z * z * R_Calculator.Ri(Temp, y) / satz / EMDEvo.Get_Hubble_at_T_Period(EMD::EPE, Temp1) - y * entropy_term;
    return res;
}

VD LeptogenesisBE::dYdX3(REAL x, VD y)
{
    REAL z = exp(x);
    VD res(DOF);
    REAL Temp = R_Calculator.Get_Temperature_from_Z(z);
    REAL Temp1 = R_Calculator.Get_Temperature_from_Z(1);
    res = z * z * R_Calculator.Ri(Temp, y) / Entropy_Density(Temp) / EMDEvo.Get_Hubble_at_T_Period(EMD::RDE, Temp1);
    return res;
}

void LeptogenesisBE::Solve()
{
    R_Calculator.Update();
    // * Starting solving the Boltzmann Equation;
    zi[0] = R_Calculator.Get_Z_from_Temperature(EMDEvo.Get_TInflation());
    zi[1] = R_Calculator.Get_Z_from_Temperature(EMDEvo.Get_Ti());
    zi[2] = R_Calculator.Get_Z_from_Temperature(EMDEvo.Get_Te());
    zi[3] = R_Calculator.Get_Z_from_Temperature(EMDEvo.Get_Tr());
    zi[4] = R_Calculator.Get_Z_from_Temperature(EMDEvo.Get_Tf());

    double TSphaleron = 100;
    double zsphaleron = R_Calculator.Get_Z_from_Temperature(TSphaleron);
    VD Y_PERIOD_BEGIN(DOF);
    for (int i = 0; i < DOF; i++)
    {
        Y_PERIOD_BEGIN[i] = 0;
    }
    for (int i = 0; i < NPeriods; i++)
    {
        PeriodID = i;
        Set_X_BEGIN(log(zi[i]));
        if (zi[i + 1] < zsphaleron)
        {
            Set_X_END(log(zi[i + 1]));
            Set_BOUNDARY_CONDITION(Y_PERIOD_BEGIN);
            solver.Solve(0.01);
            logz.push_back(solver.Get_Solution_X());
            Yields.push_back(solver.Get_Solution_Y());
            Y_PERIOD_BEGIN = solver.Get_Solution_Y_END();
        }
        else
        {
            Set_X_END(log(zsphaleron));
            Set_BOUNDARY_CONDITION(Y_PERIOD_BEGIN);
            solver.Solve(0.01);
            logz.push_back(solver.Get_Solution_X());
            Yields.push_back(solver.Get_Solution_Y());
            Y_PERIOD_BEGIN = solver.Get_Solution_Y_END();
            Set_X_BEGIN(log(zsphaleron));
            Set_X_END(log(zi[i + 1]));
            Y_PERIOD_BEGIN[0] = 0;
            Set_BOUNDARY_CONDITION(Y_PERIOD_BEGIN);
            solver.Solve(0.01);
            logz.push_back(solver.Get_Solution_X());
            Yields.push_back(solver.Get_Solution_Y());
            Y_PERIOD_BEGIN = solver.Get_Solution_Y_END();
        }
    }
}

void LeptogenesisBE::Dump_Solution(string filename)
{
    ofstream output(filename.c_str());
    output << "z\t";
    for (size_t i = 0; i < DOF; i++)
    {
        output << "y_" << i << "\t";
    }
    output << endl;
    output << scientific;
    output << showpos;
    output << setprecision(8);
    for (int i = 0; i < logz.size(); i++)
    {
        for (int j = 0; j < logz[i].size(); j++)
        {
            output << exp(logz[i][j]) << "\t";
            for (int k = 0; k < DOF; k++)
            {
                output << Yields[i][j][k] << "\t";
            }
            output << endl;
        }
    }
    output.close();
}
