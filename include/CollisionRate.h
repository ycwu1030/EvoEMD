#ifndef _COLLISION_RATE_H_
#define _COLLISION_RATE_H_

#include "Particles.h"
#include "RealTypes.h"

struct Process_Amp {
    // * We consider 2 -> 2 processes or 1 -> 2 decays
    // * 0 1 2 for s-channel, t-channel, u-channel
    // * For one single diagram, it has three cases: s, t, u channel (contact is counted as s channel)
    // * For squared, then it has 3x3 cases.
    // * But in real calculation, we will add s-t and t-s together, so Amp_xxx[0][1] == s-t + t-s = Amp_xxx[1][0] etc
    // * For decay process, treated as s-channel, [0][0] is used.
    // * The last index is for store expression in cth expansion, [x][y][0] + [x][y][1]*cth + [x][y][2]*cth^2
    REAL Amp_Numerator[3][3][3];
    REAL Amp_Denominator[3][3][3];
};

class Amplitude {
    // * Amplitude used to calculate the collision rate;
    // * As the purpose of this class is to assist the collision rate calculation and solving Boltzmann equation
    // * The relevant dof is summed insided Get_Amp()
    // * The dof has two parts:
    // * 1. The spin dof: the amp is summed over spin dof for all external particles
    // * 2. The species dof: CP conjugated component is summed as 1 -> 2 3 is added with 1bar -> 2bar 3bar etc.
    // *                     If CP violating rate is needed, (1->23) - (1bar->2bar3bar) is provided
public:
    typedef std::vector<Pseudo_Particle *> INITIAL_STATES;
    typedef std::vector<Pseudo_Particle *> FINAL_STATES;

    Amplitude();
    virtual ~Amplitude();
    virtual Process_Amp Get_Amp() = 0;
    virtual REAL Get_Coeff() = 0;

    unsigned N_INITIAL;
    INITIAL_STATES INITIAL;

    unsigned N_FINAL;
    FINAL_STATES FINAL;
};

class Collision_Rate {
    // For any process x -> y
    // CP conserving rate means, gamma(x->y) + gamma(xbar -> ybar)
    // CP violating rate means, gamma(x->y) - gamma(xbar -> ybar)
protected:
    // * This extra factor should be added for particular process,
    // * if one omit the internal symmetry for simplicity during the calculation
    // * e.g. In SM, L and Phi are SU2 doublet
    // * In some case, we don't need to manipulate the component
    // * but just treat them as a whole object.
    // * In this case, extra factor should be added.
    double Extra_Factor_for_Internal_Symmetry;

    // * ptr to Amplitude class, but we don't own it.
    Amplitude *amp;

public:
    Collision_Rate(Amplitude *amp) { this->amp = amp; }
    virtual ~Collision_Rate(){};

    virtual REAL Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat) = 0;
    virtual REAL Get_Collision_Rate(REAL T) = 0;
};

class Decay12_Rate : public Collision_Rate {
    // * For two body decay
public:
    Decay12_Rate(Amplitude *amp) : Collision_Rate(amp){};
    ~Decay12_Rate(){};

    virtual REAL Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat);
    virtual REAL Get_Collision_Rate(REAL T);
};

class Scatter22_Rate : public Collision_Rate {
protected:
public:
    Scatter22_Rate(Amplitude *amp) : Collision_Rate(amp){};
    ~Scatter22_Rate(){};

    virtual REAL Get_Collision_Rate(REAL T);
};

class Process {
public:
    using INITIAL_STATES = Amplitude::INITIAL_STATES;
    using FINAL_STATES = Amplitude::FINAL_STATES;

    Process(Amplitude *amp);
    ~Process();

    REAL Get_CP_Conserving_Rate(REAL T);
    REAL Get_CP_Conserving_Coeff(REAL T);
    REAL Get_CP_Violating_Rate(REAL T);
    REAL Get_CP_Violating_Coeff(REAL T);

protected:
    INITIAL_STATES INIT;
    FINAL_STATES FINAL;
    Amplitude *amp;
    Collision_Rate *CR_Calculator;
};

#endif  //_COLLISION_RATE_H_
