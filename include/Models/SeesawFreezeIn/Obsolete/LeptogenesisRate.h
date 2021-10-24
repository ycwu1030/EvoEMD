#ifndef __LEPTOGENESIS_RATE_H__
#define __LEPTOGENESIS_RATE_H__

#include "EvoEMD/EMD.h"
#include "EvoEMD/Neutrino.h"

// int gamma_Integrand(const int *ndim, const REAL x[], const int *ncomp,
//                     REAL ff[], void *_params);

namespace EvoEMD {

// * Calculation Collision Rate for Boltzmann
// * We are considering 5 components:
// * 0. Lepton Number
// * 1. N1, the first heavy neutrino
// * 2. N2, the second heavy neutrino
// * 3. Chi, the fermion in dark sector
// * 4. S, the scalar in dark sector
class LeptogenesisRate {
public:
    LeptogenesisRate();
    ~LeptogenesisRate(){};

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
     * @brief  Update parameters that do not depend on temperature
     * @note   Should be called prior solving corresponding Boltzmann equation
     * @retval None
     */
    void Update();

    /**
     * @brief  Get the DOF of current system
     * @note   It is 5.
     * @retval
     */
    int Get_DOF() { return DOF; }

    /**
     * @brief  Calculate the collision rate for Boltzmann equation;
     * @note
     * @param  Temp: The temperature in GeV
     * @param  y: The yields for the 5 components
     * @retval
     */
    VD Ri(REAL Temp, VD y);

    /**
     * @brief  Convert from z = M/T to T
     * @note
     * @param  z: z = M/T
     * @retval
     */
    REAL Get_Temperature_from_Z(REAL z) { return MNR1 / z; }
    REAL Get_Z_from_Temperature(REAL T) { return MNR1 / T; }

private:
    int DOF;
    VD Yeq;                 // * Yield at equilibrium
    REAL MNR1, MNR2, MNR3;  // * Heavy Neutrino masses
    REAL MCHI, MS;          // * Dark Sector masses
    REAL LamX;              // * The couplings in dark sector

    REAL GammaN1, GammaN2;  // * The width for N1 and N2
    REAL eps1, eps2;        // * The CP-asymmetry in N1 and N2 decays

    Nu_TypeI_SeeSaw Nu_Param;  // * The Type-I seesaw neutrino system

    // * The gamma's for several processes
    // * The index is for Heavy Neutrino involved in the process
    REAL gamma_Ni_Chi_S[2];
    REAL gamma_S_Chi_Ni[2];
    REAL gamma_Chi_S_Ni[2];
    REAL gamma_L_Phi_Chi_S;
    REAL gamma_Ni_Ni_Chi_Chi[2][2];
    REAL gamma_Ni_Ni_S_S[2][2];
    REAL gamma_Ni_L_Phi[2];

    /**
     * @brief  Calculate the width for N1 and N1
     * @note   Will be called in Update();
     * @retval None
     */
    void Calc_N_Width();

    /**
     * @brief Calculate the CP asymmetry for N1 and N2
     * @note
     * @retval None
     */
    void Calc_CP_Asym();

    /**
     * @brief  Calling following individual function to calculate gammas for different processes.
     * @note
     * @param  Temp: Temperature
     * @retval None
     */
    void Calc_Gammas(REAL Temp);

    /**
     * @brief  Calculate gamma for N_i->Chi S decay.
     * @note
     * @param  Temp: Temperature
     * @param  i: The index for heavy neutrino
     * @retval
     */
    REAL Calc_NChiS_Gamma(REAL Temp, int i);

    /**
     * @brief  Calculate gamma for Chi -> S N_i decay
     * @note
     * @param  Temp: Temperature
     * @param  i: The index for heavy neutrino
     * @retval
     */
    REAL Calc_ChiSN_Gamma(REAL Temp, int i);

    /**
     * @brief  Calculate gamma for S -> Chi N_i decay
     * @note
     * @param  Temp: Temperature
     * @param  i: The index for heavy neutrino
     * @retval
     */
    REAL Calc_SChiN_Gamma(REAL Temp, int i);

    /**
     * @brief  Calculate gamma for N_i -> L Phi decay
     * @note
     * @param  Temp: Temperature
     * @param  i: The index for heavy neutrino
     * @retval
     */
    REAL Calc_NLPhi_Gamma(REAL Temp, int i);

    /**
     * @brief  Calculate gamma for L Phi -> Chi S
     * @note
     * @param  Temp: Temperature
     * @retval
     */
    REAL Calc_LPhiChiS_Gamma(REAL Temp);

    /**
     * @brief  Calculate gamma for L Phi -> Lbar Phibar
     * @note
     * @param  Temp: Temperature
     * @retval
     */
    REAL Calc_LPhiLbarPhibar_Gamma(REAL Temp);

    /**
     * @brief  Calculate gamma for N_i N_j -> Chi Chi
     * @note
     * @param  Temp: Temperature
     * @param  i: The index for heavy neutrino
     * @param  j: The index for the second heavy neutrino
     * @retval
     */
    REAL Calc_NNChiChi_Gamma(REAL Temp, int i, int j);

    /**
     * @brief  Calculate gamma for N_i N_j -> S S
     * @note
     * @param  Temp: Temperature
     * @param  i: The index for heavy neutrino
     * @param  j: The index for the second heavy neutrino
     * @retval
     */
    REAL Calc_NNSS_Gamma(REAL Temp, int i, int j);

    /**
     * @brief  Calculate |M|^2 integrate over solid angle of the final states times Kallen factors for both initial and
     * final states
     * @note
     * @param  s: The C.O.M energy square
     * @param  processid: process id, 0 for L Phi -> Chi S, 1 for N N -> Chi Chi, 2 for N N -> S S
     * @param  i: The index for heavy neutrino (no use for processid=0)
     * @param  j: The index for second heavy neutrino (no use for processid=0)
     * @retval
     */
    REAL SqAmp_dOmega_with_Kallen(REAL s, int processid, int i, int j);
    REAL SqAmp_dOmega_with_Kallen_LPhiChiS(REAL s);
    REAL SqAmp_dOmega_with_Kallen_NNChiChi(REAL s, int i, int j);
    REAL SqAmp_dOmega_with_Kallen_NNSS(REAL s, int i, int j);

    REAL SqAmp_dOmega_with_Kallen_NNChiChi_v2(REAL s, int i, int j);
    REAL SqAmp_dOmega_with_Kallen_NNSS_v2(REAL s, int i, int j);
    REAL SqAmp_dOmega_with_Kallen_LPhiLbarPhiBar(REAL s);

    /**
     * @brief  Calculate the Yield at equilibrium
     * @note
     * @param  Temp: Temperature
     * @retval None
     */
    void Calc_Yield_Equilibrium(REAL Temp);

    REAL R0(REAL Temp, VD y);  // For Lepton Number;
    REAL R1(REAL Temp, VD y);  // For N1
    REAL R2(REAL Temp, VD y);  // For N2
    REAL R3(REAL Temp, VD y);  // For Chi
    REAL R4(REAL Temp, VD y);  // For S

    // * This is the function that feed to Cuba for integral.
    // * This function did nothing but calling the corresponding SqAmp_dOmega_with_Kallen function
    friend int gamma_Integrand(const int *ndim, const REAL x[], const int *ncomp, REAL ff[], void *_params);
};
}  // namespace EvoEMD
#endif  //__LEPTOGENESIS_RATE_H__
