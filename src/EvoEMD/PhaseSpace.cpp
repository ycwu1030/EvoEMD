#include "EvoEMD/PhaseSpace.h"

#include <cmath>

namespace EvoEMD {
REAL Kallen_Lam(REAL x, REAL y, REAL z) {
    REAL lam = x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x;
    return lam;
}

REAL Ei(REAL sqrt_s, REAL mi, REAL mj) {
    REAL s = sqrt_s * sqrt_s;
    REAL mi2 = mi * mi;
    REAL mj2 = mj * mj;
    REAL Energy = (s + mi2 - mj2) / 2 / sqrt_s;
    return Energy;
}

REAL Pi(REAL sqrt_s, REAL mi, REAL mj) {
    REAL s = sqrt_s * sqrt_s;
    REAL mi2 = mi * mi;
    REAL mj2 = mj * mj;
    REAL lam = Kallen_Lam(1.0, mi2 / s, mj2 / s);
    REAL Momentum = sqrt_s / 2 * sqrt(lam);
    return Momentum;
}
}  // namespace EvoEMD
