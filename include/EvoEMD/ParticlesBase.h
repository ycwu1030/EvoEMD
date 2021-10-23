#ifndef _PARTICLE_BASE_H_
#define _PARTICLE_BASE_H_

#include <set>

#include "EvoEMD/PhysicsParameters.h"
#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

class Process;
class Particle_Client {
public:
    Particle_Client(){};
    virtual ~Particle_Client(){};

    virtual void Update_Particle_Info() = 0;
};

class Particle_Base {
protected:
    // * Objects that use particle information
    // * In any time, particle information is updated, the client should update their own properties
    std::set<Particle_Client *> Particle_Client_Set;

    // * Processes that involve current particle
    // * The process calculates its value lazily, it will acquire particle's infor when it needs
    // * This set is used to keep aware what are the processes involving current particle
    std::set<Process *> Process_Set;

public:
    Particle_Base(){};
    virtual ~Particle_Base(){};

    void Register_Client(Particle_Client *);
    void Register_Process(Process *);

    void Notify_Client();
};

class Pseudo_Particle : public Particle_Base {
    // * Pseudo-Particle: not really a particle
    // * Stored Pseudo-Particle information:
    // * Masses, whether Massless
    // * PID: particle id;
    // * DOF: degree of freedom, (particle and antiparticle count seperately)
    // * Calculate Number Density or Yield at Equilibrium;
protected:
    bool selfconjugate;
    bool massless;
    bool thermalized;
    Base_Parameter *p_mass;
    Base_Parameter *p_width;
    int PID;
    int DOF;
    REAL Get_Equilibrium_Number_Density_per_DOF_Maxwell(const REAL T) const;
    REAL Get_Equilibrium_Yield_per_DOF_Maxwell(const REAL T) const;
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const = 0;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const = 0;

public:
    /**
     * @brief
     * @note
     * @param  PID: The PID for a particle
     * @param  DOF: The degree of freedom of the particle, particle and antiparticle will be counted seperately
     *              On the other hand, DOF is the one used to calculate the equilibrium number density which will be
     *              used in Boltzmann equation. So be sure this DOF is consistent with the collision rate.
     * @param  mass: Pointer to mass parameter, if nullptr (default), it is assumed the particle is massless
     * @param  width: Pointer to width parameter, if nullptr (default), it is assumed the particle is stable
     * @param  selfconjugate: Whether
     * @retval
     */
    Pseudo_Particle(int PID, int DOF, Base_Parameter *mass = nullptr, Base_Parameter *width = nullptr,
                    bool selfconjugate = false);
    virtual ~Pseudo_Particle(){};

    int Get_PID() const { return PID; }
    int Get_DOF() const { return DOF; }
    double Get_Mass() const { return p_mass->Get_Value(); }
    bool Is_Massless() const { return massless; }
    bool Is_Selfconjugate() const { return selfconjugate; }

    REAL Get_Equilibrium_Number_Density_at_T(const REAL T) const {
        return DOF * Get_Equilibrium_Number_Density_per_DOF(T);
    }
    REAL Get_Equilibrium_Yield_at_T(const REAL T) const { return DOF * Get_Equilibrium_Yield_per_DOF(T); };

    REAL Yield;
    REAL Numer_Density;
    void Set_Mass(double mass);
};

class Fermion : public Pseudo_Particle {
protected:
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const override;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const override;

public:
    Fermion(int PID, int DOF, Base_Parameter *mass = nullptr, Base_Parameter *width = nullptr,
            bool selfconjugate = false);
    ~Fermion(){};
};

class Boson : public Pseudo_Particle {
protected:
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const override;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const override;

public:
    Boson(int PID, int DOF, Base_Parameter *mass = nullptr, Base_Parameter *width = nullptr,
          bool selfconjugate = false);
    ~Boson(){};
};

}  // namespace EvoEMD
#endif  //_PARTICLE_H_
