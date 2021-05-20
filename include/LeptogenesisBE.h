#ifndef __LEPTOGENESIS_BEQ_H__
#define __LEPTOGENESIS_BEQ_H__

#include "RungeKutta.h"
#include "EMD.h"
#include "Neutrino.h"
#include "LeptogenesisRate.h"

// Y0 for Lepton Number, Y1 for N1, Y2 for N2, Y3 for Chi, Y4 for S
class LeptogenesisBE : public ODE_FUNCS
{
public:
    LeptogenesisBE();
    ~LeptogenesisBE(){};

    void Set_Temperatures(REAL Ti, REAL Tr, REAL Tf = 1e-3, REAL TInflation = 1e15);
    void Set_Nu_Masses(REAL M1, REAL M2, REAL M3);
    void Set_Nu_MassOrdering(Nu_TypeI_SeeSaw::MassOrdering od = Nu_TypeI_SeeSaw::NORMAL_ORDER);
    void Set_Dark_Sector_Masses(REAL mchi, REAL ms);
    void Set_LambdaX(REAL lamx);
    virtual VD dYdX(REAL x, VD y);

    void Solve();
    void Dump_Solution(std::string filename);

private:
    static const int NPeriods = 4;
    static const int NPoints = NPeriods + 1;

    int PeriodID;

    EMD EMDEvo;
    LeptogenesisRate R_Calculator;
    RungeKutta solver;

    REAL zi[NPoints];

    VD dYdX0(REAL x, VD y); // For ERD;
    VD dYdX1(REAL x, VD y); // For EMD;
    VD dYdX2(REAL x, VD y); // For EP;
    VD dYdX3(REAL x, VD y); // For RD;

    std::vector<VD> logz;
    std::vector<VVD> Yields;
};

#endif //__LEPTOGENESIS_BEQ_H__
