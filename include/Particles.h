#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <set>

#include "RealTypes.h"

class Process;
class Particle {
protected:
    bool selfconjugate;
    bool massless;
    double mass;
    int PID;
    int DOF;
    REAL Get_Equilibrium_Number_Density_per_DOF_Maxwell(const REAL T) const;
    REAL Get_Equilibrium_Yield_per_DOF_Maxwell(const REAL T) const;
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const = 0;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const = 0;
    std::set<Process *> Process_List;

public:
    /**
     * @brief  Constructor for massless particle
     * @note
     * @param  PID: particle ID, unique and positive, the corresponding negative one is for the antiparticle
     * @param  DOF: The internal dof of the particle, only for particle (antiparticle is not included)
     * @param  selfconjugate: Self conjugate or not
     * @retval
     */
    Particle(int PID, int DOF, bool selfconjugate = false);

    /**
     * @brief  Constructor for massive particle
     * @note
     * @param  mass: mass for the particle
     * @param  PID: particle ID, unique and positive, the corresponding negative one is for the antiparticle
     * @param  DOF: The internal dof of the particle, only for particle (antiparticle is not included)
     * @param  selfconjugate: Self conjugate or not
     * @retval
     */
    Particle(double mass, int PID, int DOF, bool selfconjugate = false);
    virtual ~Particle(){};

    int Get_PID() const { return PID; }
    double Get_Mass() const { return mass; }
    bool Is_Massless() const { return massless; }
    bool Is_Selfconjugate() const { return selfconjugate; }
    virtual int Get_DOF() const { return DOF; }
    REAL Get_Equilibrium_Number_Density_at_T(const REAL T) const {
        return DOF * Get_Equilibrium_Number_Density_per_DOF(T);
    }
    REAL Get_Equilibrium_Yield_at_T(const REAL T) const { return DOF * Get_Equilibrium_Yield_per_DOF(T); };

    void Register_Process(Process *proc) { Process_List.insert(proc); }
    void Set_Mass(double mass);
};

class Fermion : public Particle {
protected:
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const override;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const override;

public:
    Fermion(int PID, int DOF, bool selfconjugate = false);
    Fermion(double mass, int PID, int DOF, bool selfconjugate = false);
    ~Fermion(){};
};

class Boson : public Particle {
protected:
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const override;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const override;

public:
    Boson(int PID, int DOF, bool selfconjugate = false);
    Boson(double mass, int PID, int DOF, bool selfconjugate = false);
    ~Boson(){};
};

#endif  //_PARTICLE_H_
