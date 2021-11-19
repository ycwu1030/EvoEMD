#ifndef _COLLISION_RATE_H_
#define _COLLISION_RATE_H_

#include <map>

#include "EvoEMD/ParticleBase.h"
#include "EvoEMD/RealTypes.h"

namespace EvoEMD {
// struct Process_Amp_Single_Diagram {
//     typedef VD Result_Type;
//     typedef Result_Type Numerator_Type;
//     typedef int Propagator_ID;
//     typedef std::pair<Propagator_ID, Result_Type> Propagator_Type;
//     typedef std::pair<Propagator_Type, Propagator_Type> Denominator_Type;
//     Numerator_Type Numerator;
//     Denominator_Type Denominator;
// };
// typedef std::vector<Process_Amp_Single_Diagram> Process_Amp;

// template <typename... T>
// inline VD Build_Numerator(T... args) {
// return VD{args...};
// }
// inline VD Build_Numerator(REAL n0, REAL n1, REAL n2) { return VD{n0, n1, n2}; }

// template <typename... T>
// inline Process_Amp_Single_Diagram::Propagator_Type Build_Propagator(int ID, T... args) {
//     return std::make_pair(ID, VD{args...});
// }
// inline Process_Amp_Single_Diagram::Propagator_Type Build_Propagator(int ID, REAL d0, REAL d1) {
//     return d1 == 0 ? std::make_pair(ID, VD{d0}) : std::make_pair(ID, VD{d0, d1});
// }
// inline Process_Amp_Single_Diagram::Propagator_Type Build_Propagator(int ID, REAL d0) {
//     return std::make_pair(ID, VD{d0});
// }

// inline Process_Amp_Single_Diagram::Denominator_Type Build_Denominator(Process_Amp_Single_Diagram::Propagator_Type p1,
//                                                                       Process_Amp_Single_Diagram::Propagator_Type p2)
//                                                                       {
//     return std::make_pair(p1, p2);
// }

class Amplitude_Base : public Parameter_Base {
    // * Amplitude_Base used to calculate the collision rate;
    // * As the purpose of this class is to assist the collision rate calculation and solving Boltzmann equation
    // * The relevant dof is summed insided Get_Amp()
    // * The dof has two parts:
    // * 1. The spin dof: the amp is summed over spin dof for all external particles
    // * 2. The species dof: CP conjugated component is summed as 1 -> 2 3 is added with 1bar -> 2bar 3bar etc.
    // *                     If CP violating rate is needed, (1->23) - (1bar->2bar3bar) is provided
    // *                     Flavor dof is also summed
    // * Amplitude_Base is also derived from Parameter_Base, such that Amplitude_Base will update any sqrt(s)
    // * independent parameter inside (In Update_Value(input)).
    // * Update_Amp(sqrt_shat) will update the amplitdue.
    // * Currently only 1->2 and 2->2 cases.
    // * Get_Amp(REAL sqrt_shat) return the \int d\Omega |M|^2
public:
    typedef std::vector<Particle_Base *> INITIAL_STATES;
    typedef std::vector<Particle_Base *> FINAL_STATES;

    Amplitude_Base(std::string name) : Parameter_Base(name) {}
    virtual ~Amplitude_Base(){};
    virtual void Update_Value(REAL input) = 0;
    virtual void Update_Amp(REAL sqrt_shat) = 0;
    REAL Get_Amp(REAL sqrt_shat) {
        Update_Amp(sqrt_shat);
        return amp_res;
    };
    virtual REAL Get_Coeff(REAL T, int PID) = 0;

    unsigned N_INITIAL;
    INITIAL_STATES INITIAL;

    unsigned N_FINAL;
    FINAL_STATES FINAL;

protected:
    REAL amp_res;
};

class Collision_Rate {
    // For any process x -> y
    // CP conserving rate means, gamma(x->y) + gamma(xbar -> ybar)
    // CP violating rate means, gamma(x->y) - gamma(xbar -> ybar)
protected:
    // * ptr to Amplitude_Base class, but we don't own it.
    Amplitude_Base *amp;

public:
    Collision_Rate(Amplitude_Base *amp) { this->amp = amp; }
    virtual ~Collision_Rate(){};

    // virtual REAL Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) = 0;
    REAL Get_Amp_With_Solid_Angle(REAL sqrt_shat) { return amp->Get_Amp(sqrt_shat); }
    virtual REAL Get_Collision_Rate(REAL T) = 0;
};

class Decay12_Rate : public Collision_Rate {
    // * For two body decay
public:
    Decay12_Rate(Amplitude_Base *amp) : Collision_Rate(amp){};
    ~Decay12_Rate(){};

    // virtual REAL Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat);
    virtual REAL Get_Collision_Rate(REAL T);
};

class Scatter22_Rate : public Collision_Rate {
    // protected:
    //     REAL Get_Amp_Integrate_over_Phase_Space_Single_Channel(const Process_Amp_Single_Diagram &amp_single);

public:
    Scatter22_Rate(Amplitude_Base *amp) : Collision_Rate(amp){};
    ~Scatter22_Rate(){};

    // virtual REAL Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat);
    virtual REAL Get_Collision_Rate(REAL T);
};

}  // namespace EvoEMD
#endif  //_COLLISION_RATE_H_
