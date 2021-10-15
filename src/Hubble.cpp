#include "Hubble.h"

#include <cmath>

#include "EffDOF.h"
#include "Physics_Constants.h"

Hubble_For_Single_Period::Hubble_For_Single_Period(double Ti, double Tf, bool Isentropic_in, double beta_T_in,
                                                   double beta_s_in)
    : T_start(Ti), T_end(Tf), Isentropic(Isentropic_in), beta_T(beta_T_in), beta_s(beta_s_in) {}

REAL Hubble_For_Single_Period::Get_Hubble_For_RD(REAL T) {
    REAL geT = ge(T);
    return M_PI / 3.0 * sqrt(geT / 10.0) * T * T / PHY_MP / PHY_MP;
}

Hubble_RD::Hubble_RD(REAL T_start, REAL T_end) : Hubble_For_Single_Period(T_start, T_end) {}

REAL Hubble_RD::Get_Hubble_at_T(REAL T) { return Get_Hubble_For_RD(T); }

Hubble_EMD::Hubble_EMD(REAL T_start, REAL T_end) : Hubble_For_Single_Period(T_start, T_end) {
    HRD_at_T_start = Get_Hubble_For_RD(T_start);
}

REAL Hubble_EMD::Get_Hubble_at_T(REAL T) { return HRD_at_T_start * pow(T / T_start, 1.5); }

Hubble_EP::Hubble_EP(REAL T_start, REAL T_end) : Hubble_For_Single_Period(T_start, T_end, false, 3.0 / 8.0, 9.0 / 8.0) {
    HRD_at_T_end = Get_Hubble_For_RD(T_end);
    ge_at_T_end = ge(T_end);
}

REAL Hubble_EP::Get_Hubble_at_T(REAL T) {
    REAL geT = ge(T);
    return HRD_at_T_end * geT / ge_at_T_end * pow(T / T_end, 4);
}
