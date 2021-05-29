#ifndef __LEPTOGENESIS_BEQ_H__
#define __LEPTOGENESIS_BEQ_H__

#include "EMD.h"
#include "LeptogenesisRate.h"
#include "Neutrino.h"
#include "RungeKutta.h"

// * Y0 for Lepton Number, Y1 for N1, Y2 for N2, Y3 for Chi, Y4 for S
// * This class assemble useful object to build up the Boltzmann equation
class LeptogenesisBE : public ODE_FUNCS {
public:
    LeptogenesisBE();
    ~LeptogenesisBE(){};

    /**
     * @brief  Set up the temperature for the Cosmology History
     * @note   This function will call EMD.Set_Temperature
     * @param  Ti: Initial temperature of EMD
     * @param  Tr: The end temperature of EMD
     * @param  Tf: The lowest temperature we will consider
     * @param  TInflation: The highest temperature we will consider
     * @retval None
     */
    void Set_Temperatures(REAL Ti, REAL Tr, REAL Tf = 1e-3, REAL TInflation = 1e15);

    /**
     * @brief  Set the Heavy Neutrino masses
     * @note   It is better to have M1<=M2<=M3
     *         This will call LeptogenesisRate functions
     * @param  M1: Mass of one NR
     * @param  M2: Mass of the second NR
     * @param  M3: Mass of the last NR
     * @retval None
     */
    void Set_Nu_Masses(REAL M1, REAL M2, REAL M3);

    /**
     * @brief  Set the light neutrino mass ordering
     * @note
     * @param  od: NORMAL_ORDER or INVERTED_ORDER
     * @retval None
     */
    void Set_Nu_MassOrdering(Nu_TypeI_SeeSaw::MassOrdering od = Nu_TypeI_SeeSaw::NORMAL_ORDER);

    /**
     * @brief  Set the dark sector masses
     * @note
     * @param  mchi: mass for chi
     * @param  ms: mass for s
     * @retval None
     */
    void Set_Dark_Sector_Masses(REAL mchi, REAL ms);

    /**
     * @brief  Set the dark sector couplings
     * @note
     * @param  lamx: the coupling
     * @retval None
     */
    void Set_LambdaX(REAL lamx);

    /**
     * @brief  Provide the dY/dX for the RungeKutta Solver
     * @note
     * @param  x:
     * @param  y:
     * @retval RETURN dy/dx at x
     */
    virtual VD dYdX(REAL x, VD y) override;

    /**
     * @brief  Calling to solve the Boltzmann Equation
     * @note
     * @retval None
     */
    void Solve();

    /**
     * @brief  Dump the solution into a file
     * @note
     * @param  filename: the file name
     * @retval None
     */
    void Dump_Solution(std::string filename);

private:
    static const int NPeriods = 4;            // * We have 4 periods for this calculations ERD, EMD, EP, RD
    static const int NPoints = NPeriods + 1;  // * With four period, we will have five temperature or z
                                              // * = m/T to determine the interval of these periods

    int PeriodID;

    EMD EMDEvo;                     // * Using EMD to get Hubble parameters
    LeptogenesisRate R_Calculator;  // * Calculate the rate for Boltzmann equation
    RungeKutta solver;              // * RungeKutta solver to get the evolutions

    REAL zi[NPoints];  // * Store the end points determining the four periods.

    // * dYdX for four different periods.
    VD dYdX0(REAL x, VD y);  // For ERD;
    VD dYdX1(REAL x, VD y);  // For EMD;
    VD dYdX2(REAL x, VD y);  // For EP;
    VD dYdX3(REAL x, VD y);  // For RD;

    // * vectors storing the results of the Boltzmann equation
    std::vector<VD> logz;
    std::vector<VVD> Yields;
};

#endif  //__LEPTOGENESIS_BEQ_H__
