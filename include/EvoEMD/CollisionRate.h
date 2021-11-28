#ifndef _COLLISION_RATE_H_
#define _COLLISION_RATE_H_

#include <map>

#include "EvoEMD/ParticleBase.h"
#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

class Amplitude_Base : public Parameter_Base {
    // * Amplitude_Base used to calculate the collision rate;
    // * As the purpose of this class is to assist the collision rate calculation and solving Boltzmann equation
    // * The relevant dof is summed insided Get_Amp()
    // * The dof has two parts:
    // * 1. The spin dof: the amp is summed over spin dof for all external particles
    // * 2. The species dof: CP conjugated component is summed as 1 -> 2 3 is added with 1bar -> 2bar 3bar etc.
    // *                     If CP violating rate is needed, (1->23) - (1bar->2bar3bar) is provided
    // *                     Flavor dof is also summed
    // * Amplitude_Base is also inherited from Parameter_Base, such that Amplitude_Base will update any sqrt(s)
    // * independent parameter inside (In Update_Value(input)).
    // * Update_Amp(sqrt_shat) will update the amplitdue.
    // * Get_Amp(REAL sqrt_shat) returns |M|^2 integrated over final state phase space
public:
    typedef std::vector<Particle_Base *> INITIAL_STATES;
    typedef std::vector<Particle_Base *> FINAL_STATES;

    Amplitude_Base(std::string name) : Parameter_Base(name) {}
    virtual ~Amplitude_Base(){};

    /**
     * @brief Update any sqrt_shat independent values that will be used to compute amplitude
     * @note
     * @param  input: dummy input
     * @retval None
     */
    virtual void Update_Value(REAL input){};

    /**
     * @brief  Update the offset used to build BE
     * @note
     * @param  T: the temperature at which to evaluate the coeff
     * @param  PID: the PID for the particle we want to evaluate the coeff
     * @retval
     */
    virtual REAL Get_Offset(REAL T, int PID) { return 0; }

    /**
     * @brief  Update the amplitude (integrated over final state phase space), and store the result in amp_res
     * @note
     * @param  sqrt_shat: the process energy
     * @retval
     */
    virtual void Update_Amp(REAL sqrt_shat) { amp_res = 0; }

    /**
     * @brief  Get the amplitude integrated over final state phase space
     * @note   Inside this function, the internal parameters will be updated first
     * @param  sqrt_shat: the process energy
     * @retval Return the amplitude integrated over final state phase space
     */
    REAL Get_Amp(REAL sqrt_shat) {
        Get_Value();
        Update_Amp(sqrt_shat);
        return amp_res;
    }

    unsigned N_INITIAL;
    INITIAL_STATES INITIAL;

    unsigned N_FINAL;
    FINAL_STATES FINAL;

protected:
    REAL amp_res;
};

class Collision_Rate {
    // For any process x -> y
protected:
    // * ptr to Amplitude_Base class, but we don't own it.
    Amplitude_Base *amp;

public:
    Collision_Rate(Amplitude_Base *amp) { this->amp = amp; }
    virtual ~Collision_Rate(){};

    REAL Get_Amp(REAL sqrt_shat) { return amp->Get_Amp(sqrt_shat); }
    virtual REAL Get_Collision_Rate(REAL T) = 0;
};

class Decay_Rate : public Collision_Rate {
public:
    Decay_Rate(Amplitude_Base *amp) : Collision_Rate(amp){};
    ~Decay_Rate(){};

    virtual REAL Get_Collision_Rate(REAL T);
};

class Scatter_Rate : public Collision_Rate {
public:
    Scatter_Rate(Amplitude_Base *amp) : Collision_Rate(amp){};
    ~Scatter_Rate(){};

    virtual REAL Get_Collision_Rate(REAL T);
};

}  // namespace EvoEMD
#endif  //_COLLISION_RATE_H_
