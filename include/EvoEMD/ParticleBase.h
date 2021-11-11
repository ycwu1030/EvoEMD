#ifndef _PARTICLE_BASE_H_
#define _PARTICLE_BASE_H_

#include <map>
#include <set>

#include "EvoEMD/ParameterBase.h"
#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

class Process;

class Particle_Base {
protected:
    // * Processes that involve current particle
    // * The process calculates its value lazily, it will acquire particle's infor when it needs
    // * This set is used to keep aware what are the processes involving current particle
    std::set<Process *> Process_Set;

public:
    Particle_Base(){};
    virtual ~Particle_Base(){};

    void Register_Process(Process *);

    std::set<Process *> Get_Process() const { return Process_Set; }
};

class Pseudo_Particle : public Particle_Base {
    // * Pseudo-Particle:
    // * Stored Pseudo-Particle information:
    // * Masses, whether Massless
    // * PID: particle id;
    // * DOF: degree of freedom, (particle and antiparticle count seperately)
    // * pseudo: whether this particle is pseudo or not
    // *    By pseudo, we mean, it is not a realistic particle, but e.g. difference between particle and anti-particle
    // *    So for pseudo, its number density/Yield can be negative.
    // * Thermalized: whether this particle is thermalized or not at the beginning of the evolution.
    // *    By default: it is assumed that all particle is in thermalization.
    // * Calculate Number Density or Yield at Equilibrium;
protected:
    const std::string name;
    bool massless;
    Parameter_Base *p_mass;
    Parameter_Base *p_width;
    const int PID;
    const int DOF;
    const bool pseudo;
    bool Thermalized;
    REAL Get_Equilibrium_Number_Density_per_DOF_Maxwell(const REAL T) const;
    REAL Get_Equilibrium_Yield_per_DOF_Maxwell(const REAL T) const;
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const = 0;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const = 0;

public:
    /**
     * @brief  ctor for a Pseudo_Particle
     * @note
     * @param  name: The name for the particle
     * @param  PID: The PID for a particle (should be unique for each particle)
     * @param  DOF: The degree of freedom of the particle, particle and antiparticle will be counted seperately
     *              On the other hand, DOF is the one used to calculate the equilibrium number density which will be
     *              used in Boltzmann equation. So be sure this DOF is consistent with the collision rate.
     * @param  mass: Pointer to mass parameter, if nullptr (default), it is assumed the particle is massless
     * @param  width: Pointer to width parameter, if nullptr (default), it is assumed the particle is stable
     * @param  pseudo: whether the particle is pseudo or not. If it is pseudo, then its equilibrium is not at Yeq, and
     * the density of it can be negative. If it is not pseudo, it is some basic particles, its equilibrium is just Yeq,
     * and the density can not be negative.
     * @retval
     */
    Pseudo_Particle(std::string name, int PID, int DOF, Parameter_Base *mass = nullptr, Parameter_Base *width = nullptr,
                    bool pseudo = false);
    virtual ~Pseudo_Particle(){};

    int Get_PID() const { return PID; }
    int Get_DOF() const { return DOF; }
    std::string Get_Name() const { return name; }
    double Get_Mass() const {
        if (massless) {
            return 0;
        } else {
            return p_mass->Get_Value();
        }
    }
    bool Is_Massless() const { return massless; }
    bool Is_Pseudo() const { return pseudo; }
    bool Get_Init_Thermal_Status() const { return Thermalized; }
    void Set_Init_Thermal_Status(bool thermal = true) { Thermalized = thermal; }

    REAL Get_Equilibrium_Number_Density_at_T(const REAL T) const {
        return DOF * Get_Equilibrium_Number_Density_per_DOF(T);
    }
    REAL Get_Equilibrium_Yield_at_T(const REAL T) const { return DOF * Get_Equilibrium_Yield_per_DOF(T); };

    REAL Yield;

    // * 1 - Yield/Yield_eq
    // * At sufficient high temperature, particle might be in thermal equilibrium
    // * In that case, this value is actually 0. But due to numerical issue, it can be a small number
    // * However, in calculating collision rate, it might be multiplied by a extremely large number,
    // * thus leads to numerical issue.
    // * So we will keep this number, especially when it is zero to avoid numerical issue.
    // * User can set their own threshold when to use this Delta_Yield_Ratio
    REAL Delta_Yield_Ratio;
    REAL Numer_Density;

    void Set_Mass(double mass);
};

class Fermion : public Pseudo_Particle {
protected:
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const override;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const override;

public:
    Fermion(std::string name, int PID, int DOF, Parameter_Base *mass = nullptr, Parameter_Base *width = nullptr,
            bool pseudo = false);
    ~Fermion(){};
};

class Boson : public Pseudo_Particle {
protected:
    virtual REAL Get_Equilibrium_Number_Density_per_DOF(const REAL T) const override;
    virtual REAL Get_Equilibrium_Yield_per_DOF(const REAL T) const override;

public:
    Boson(std::string name, int PID, int DOF, Parameter_Base *mass = nullptr, Parameter_Base *width = nullptr,
          bool pseudo = false);
    ~Boson(){};
};

class Particle_Factory {
public:
    typedef std::map<int, Pseudo_Particle *> Particle_List;

    static Particle_Factory &Get_Particle_Factory();

    Pseudo_Particle *Get_Particle(int PID);
    bool Register_Particle(Pseudo_Particle *);
    bool Register_POI(int PID);
    bool Set_Mass(int PID, double mass);
    std::set<int> Get_POI() const { return POI; }

private:
    Particle_Factory();
    ~Particle_Factory();

    Particle_List PL;   // All Particle
    std::set<int> POI;  // Particle of Interested
};

}  // namespace EvoEMD

class Register_Particle {
public:
    Register_Particle(EvoEMD::Pseudo_Particle *par) {
        EvoEMD::Particle_Factory::Get_Particle_Factory().Register_Particle(par);
    }
};

// #define REGISTER_PARTICLE(partName) Register_Particle g_register_particle_##partName(new partName)

#define REGISTER_PARTICLE(className, partName, ...)              \
    class part_##partName : public className {                   \
    public:                                                      \
        part_##partName() : className(#partName, __VA_ARGS__){}; \
    };                                                           \
    Register_Particle g_register_particle_##partName(new part_##partName)

#define RETRIVE_PARTICLE(PID) EvoEMD::Particle_Factory::Get_Particle_Factory().Get_Particle(PID)

class Register_POI {
public:
    Register_POI(int PID, bool start_with_thermal = true) {
        EvoEMD::Particle_Factory &pf = EvoEMD::Particle_Factory::Get_Particle_Factory();
        pf.Register_POI(PID);
        EvoEMD::Pseudo_Particle *pp = pf.Get_Particle(PID);
        pp->Set_Init_Thermal_Status(start_with_thermal);
    }
};

#define REGISTER_POI(PID, THERMAL) Register_POI g_register_poi_##PID(PID, THERMAL)

#endif  //_PARTICLE_H_
