#include "EvoEMD/ParticleBase.h"

#include <cmath>

#include "EvoEMD/Constants.h"
#include "EvoEMD/EffDOF.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_zeta.h"

namespace EvoEMD {
void Particle_Proc::Register_Process(Process *proc) { Process_Set.insert(proc); }

Particle_Base::Particle_Base(std::string name_in, int PID_in, int DOF_in, Parameter_Base *mass, Parameter_Base *width,
                             bool pseudo_)
    : name(name_in), DOF(DOF_in), PID(PID_in), p_mass(mass), p_width(width), pseudo(pseudo_), Thermalized(true) {
    if (!p_mass) {
        massless = true;
    } else {
        massless = false;
    }
}

REAL Particle_Base::Get_Equilibrium_Number_Density_per_DOF_Maxwell(const REAL T) const {
    REAL mass = Get_Mass();
    REAL z = mass / T;
    REAL z2K2 = z > BESSEL_Z_MAX ? 0 : z * z * gsl_sf_bessel_Kn(2, z);
    REAL neq = pow(T, 3) / 2 / M_PI / M_PI * z2K2;
    return neq;
}

REAL Particle_Base::Get_Equilibrium_Yield_per_DOF_Maxwell(const REAL T) const {
    // * We define the Yield as Y = n/s
    // * s = 2 pi^2/45 gs(T) T^3
    REAL entropy_pre_factor = 2 * M_PI * M_PI / 45 * f_gs(T);
    REAL mass = Get_Mass();
    REAL z = mass / T;
    REAL z2K2 = z > BESSEL_Z_MAX ? 0 : z * z * gsl_sf_bessel_Kn(2, z);
    REAL Yeq = z2K2 / 2 / M_PI / M_PI / entropy_pre_factor;
    return Yeq;
}

void Particle_Base::Set_Mass(double mass) {
    if (massless) {
        std::cout << "Setting mass for a massless particle, the mass is ignored" << std::endl;
    } else {
        p_mass->Set_Value(mass);
    }
}

Fermion::Fermion(std::string name, int PID, int DOF, Parameter_Base *mass, Parameter_Base *width, bool pseudo)
    : Particle_Base(name, PID, DOF, mass, width, pseudo) {}

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
        REAL entropy_pre_factor = 2 * M_PI * M_PI / 45 * f_gs(T);
        return 3.0 / 4.0 * zeta3 / M_PI / M_PI / entropy_pre_factor;
    } else {
        return Get_Equilibrium_Yield_per_DOF_Maxwell(T);
    }
}

Boson::Boson(std::string name, int PID, int DOF, Parameter_Base *mass, Parameter_Base *width, bool pseudo)
    : Particle_Base(name, PID, DOF, mass, width, pseudo) {}

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
        REAL entropy_pre_factor = 2 * M_PI * M_PI / 45 * f_gs(T);
        return zeta3 / M_PI / M_PI / entropy_pre_factor;
    } else {
        return Get_Equilibrium_Yield_per_DOF_Maxwell(T);
    }
}

Particle_Factory::Particle_Factory() {}

Particle_Factory::~Particle_Factory() {
    Particle_List::iterator iter;
    for (iter = PL.begin(); iter != PL.end(); ++iter) {
        delete iter->second;
    }
}

Particle_Factory &Particle_Factory::Get_Particle_Factory() {
    static Particle_Factory PF;
    return PF;
}

Particle_Base *Particle_Factory::Get_Particle(int PID) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        std::cout << "Cannot find particle with PID = " << PID << std::endl;
        return nullptr;
    } else {
        return iter->second;
    }
}

bool Particle_Factory::Register_Particle(Particle_Base *part) {
    if (!part) {
        return false;
    }
    int PID = part->Get_PID();
    Particle_List::iterator iter = PL.find(PID);
    if (iter != PL.end()) {
        std::cout << "There is another particle with PID = " << PID << std::endl;
        return false;
    }
    PL[PID] = part;
    return true;
}

bool Particle_Factory::Register_POI(int PID) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        std::cout << "Cannot find particle with PID = " << PID << std::endl;
        return false;
    } else {
        POI.insert(PID);
        return true;
    }
}

bool Particle_Factory::Set_Mass(int PID, double mass) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        std::cout << "Cannot find particle with PID = " << PID << std::endl;
        return false;
    } else {
        iter->second->Set_Mass(mass);
        return true;
    }
}

}  // namespace EvoEMD
