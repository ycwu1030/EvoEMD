#include "Particles.h"

#include <cmath>

#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_zeta.h"

Particle::Particle(int PID_in, int DOF_in, bool selfconjugate_in)
    : mass(0), DOF(DOF_in), PID(PID_in), selfconjugate(selfconjugate_in), massless(true) {}

Particle::Particle(double mass_in, int PID_in, int DOF_in, bool selfconjugate_in)
    : mass(mass_in), PID(PID_in), DOF(DOF_in), selfconjugate(selfconjugate_in), massless(false) {}

REAL Particle::Get_Equilibrium_Number_Density_per_DOF_Maxwell(const REAL T) const {
    REAL z = mass / T;
    REAL z2K2 = z * z * gsl_sf_bessel_Kn(2, z);
    REAL neq = pow(T, 3) / 2 / M_PI / M_PI * z2K2;
    return neq;
}

Fermion::Fermion(int PID, int DOF, bool selfconjugate) : Particle(PID, DOF, selfconjugate) {}

Fermion::Fermion(double mass, int PID, int DOF, bool selfconjugate) : Particle(mass, PID, DOF, selfconjugate) {}

REAL Fermion::Get_Equilibrium_Number_Density_per_DOF(const REAL T) const {
    static const double zeta3 = gsl_sf_zeta_int(3);
    if (massless) {
        return 3.0 * 4.0 * zeta3 / M_PI / M_PI * pow(T, 3);
    } else {
        return Get_Equilibrium_Number_Density_per_DOF_Maxwell(T);
    }
}
