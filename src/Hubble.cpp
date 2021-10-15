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
