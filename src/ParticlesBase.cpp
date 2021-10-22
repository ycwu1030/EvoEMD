#include "EvoEMD/ParticlesBase.h"

#include <cmath>

#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_zeta.h"

namespace EvoEMD {
void Particle_Base::Register_Client(Particle_Client *pc) { Particle_Client_Set.insert(pc); }
void Particle_Base::Register_Process(Process *proc) { Process_Set.insert(proc); }
void Particle_Base::Notify_Client() {
    std::set<Particle_Client *>::iterator iter = Particle_Client_Set.begin();
    for (; iter != Particle_Client_Set.end(); iter++) {
        (*iter)->Update_Particle_Info();
    }
}

Pseudo_Particle::Pseudo_Particle(int PID_in, int DOF_in, bool selfconjugate_in)
    : mass(0), DOF(DOF_in), PID(PID_in), selfconjugate(selfconjugate_in), massless(true) {}

Pseudo_Particle::Pseudo_Particle(double mass_in, int PID_in, int DOF_in, bool selfconjugate_in)
    : mass(mass_in), PID(PID_in), DOF(DOF_in), selfconjugate(selfconjugate_in), massless(false) {}

REAL Pseudo_Particle::Get_Equilibrium_Number_Density_per_DOF_Maxwell(const REAL T) const {
    REAL z = mass / T;
    REAL z2K2 = z * z * gsl_sf_bessel_Kn(2, z);
    REAL neq = pow(T, 3) / 2 / M_PI / M_PI * z2K2;
    return neq;
}

REAL Pseudo_Particle::Get_Equilibrium_Yield_per_DOF_Maxwell(const REAL T) const {
    REAL z = mass / T;
    REAL z2K2 = z * z * gsl_sf_bessel_Kn(2, z);
    REAL Yeq = z2K2 / 2 / M_PI / M_PI;
    return Yeq;
}

void Pseudo_Particle::Set_Mass(double mass) {
    if (massless) {
        std::cout << "Setting mass for a massless particle, the mass is ignored" << std::endl;
    } else {
        this->mass = mass;
        Notify_Client();
    }
}

Fermion::Fermion(int PID, int DOF, bool selfconjugate) : Pseudo_Particle(PID, DOF, selfconjugate) {}

Fermion::Fermion(double mass, int PID, int DOF, bool selfconjugate) : Pseudo_Particle(mass, PID, DOF, selfconjugate) {}

REAL Fermion::Get_Equilibrium_Number_Density_per_DOF(const REAL T) const {
    static const double zeta3 = gsl_sf_zeta_int(3);
    if (massless) {
        return 3.0 / 4.0 * zeta3 / M_PI / M_PI * pow(T, 3);
    } else {
        return Get_Equilibrium_Number_Density_per_DOF_Maxwell(T);
    }
}

REAL Fermion::Get_Equilibrium_Yield_per_DOF(const REAL T) const {
    static const double zeta3 = gsl_sf_zeta_int(3);
    if (massless) {
        return 3.0 / 4.0 * zeta3 / M_PI / M_PI;
    } else {
        return Get_Equilibrium_Yield_per_DOF_Maxwell(T);
    }
}

Boson::Boson(int PID, int DOF, bool selfconjugate) : Pseudo_Particle(PID, DOF, selfconjugate) {}

Boson::Boson(double mass, int PID, int DOF, bool selfconjugate) : Pseudo_Particle(mass, PID, DOF, selfconjugate) {}

REAL Boson::Get_Equilibrium_Number_Density_per_DOF(const REAL T) const {
    static const double zeta3 = gsl_sf_zeta_int(3);
    if (massless) {
        return zeta3 / M_PI / M_PI * pow(T, 3);
    } else {
        return Get_Equilibrium_Number_Density_per_DOF_Maxwell(T);
    }
}

REAL Boson::Get_Equilibrium_Yield_per_DOF(const REAL T) const {
    static const double zeta3 = gsl_sf_zeta_int(3);
    if (massless) {
        return zeta3 / M_PI / M_PI;
    } else {
        return Get_Equilibrium_Yield_per_DOF_Maxwell(T);
    }
}
}  // namespace EvoEMD
