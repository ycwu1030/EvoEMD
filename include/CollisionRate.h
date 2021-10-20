#ifndef _COLLISION_RATE_H_
#define _COLLISION_RATE_H_

#include <map>

#include "Particles.h"
#include "RealTypes.h"

struct Process_Amp {
    typedef enum { SS = 0, ST = 1, SU = 2, TT = 3, TU = 4, UU = 5 } Channel;
    static const Process_Amp::Channel All_Channel[6];
    typedef std::map<Process_Amp::Channel, unsigned> AMP_DIAG;
    AMP_DIAG n_diag;

    // * bool determine whether the corresponding REAL is actually exactly zero
    typedef std::pair<bool, REAL> CTH_RES;
    typedef std::vector<CTH_RES> CTH_RES_FULL;
    typedef CTH_RES_FULL NUMERATOR_STRUCTURE;
    typedef int Propagator_ID;
    typedef std::pair<Propagator_ID, CTH_RES_FULL> PROPAGATOR_STRUCTURE;
    typedef std::pair<PROPAGATOR_STRUCTURE, PROPAGATOR_STRUCTURE> DENOMINATOR_STRUCTURE;

    // * The key used to access the results: Channel + Diagram_ID;
    typedef std::pair<Channel, int> AMP_KEY;

    // * Key -> RES[3] for numerator: [0] + [1]*cth + [2]*cth^2;
    // * Key -> RES[4] for denominator: ([0] + [1]*cth)*([2]+[3]*cth)
    typedef std::map<AMP_KEY, NUMERATOR_STRUCTURE> AMP_NUM;
    typedef std::map<AMP_KEY, DENOMINATOR_STRUCTURE> AMP_DEN;
    AMP_NUM amps_numerator;
    AMP_DEN amps_denominator;

    NUMERATOR_STRUCTURE &Get_Numerator(Channel chan1 = Channel::SS, int diagram_id = 0) const;
    DENOMINATOR_STRUCTURE &Get_Denominator(Channel chan1 = Channel::SS, int diagram_id = 0) const;
};
const Process_Amp::Channel Process_Amp::All_Channel[6] = {SS, ST, SU, TT, TU, UU};

class Amplitude {
    // * Amplitude used to calculate the collision rate;
    // * As the purpose of this class is to assist the collision rate calculation and solving Boltzmann equation
    // * The relevant dof is summed insided Get_Amp()
    // * The dof has two parts:
    // * 1. The spin dof: the amp is summed over spin dof for all external particles
    // * 2. The species dof: CP conjugated component is summed as 1 -> 2 3 is added with 1bar -> 2bar 3bar etc.
    // *                     If CP violating rate is needed, (1->23) - (1bar->2bar3bar) is provided
    // *                     Flavor dof is also summed
public:
    typedef std::vector<Pseudo_Particle *> INITIAL_STATES;
    typedef std::vector<Pseudo_Particle *> FINAL_STATES;
    typedef enum { SS = 0, ST = 1, SU = 2, TT = 3, TU = 4, UU = 5 } Channel;
    typedef enum { CONST = 0, LINEAR = 1, QUAD = 2 } CTH_ORDER;

    Amplitude();
    virtual ~Amplitude(){};
    virtual void Update_Amp(REAL sqrt_shat) = 0;
    virtual const Process_Amp &Get_Amp(REAL sqrt_shat) {
        Update_Amp(sqrt_shat);
        return amp_res;
    };
    virtual REAL Get_Coeff() = 0;

    unsigned N_INITIAL;
    INITIAL_STATES INITIAL;

    unsigned N_FINAL;
    FINAL_STATES FINAL;

private:
    Process_Amp amp_res;
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
    REAL Get_Amp_Integrate_over_Phase_Space_Single_Channel(const Process_Amp::NUMERATOR_STRUCTURE &numerator,
                                                           const Process_Amp::DENOMINATOR_STRUCTURE &denominator);

public:
    Scatter22_Rate(Amplitude *amp) : Collision_Rate(amp){};
    ~Scatter22_Rate(){};

    virtual REAL Get_Amp_Integrate_over_Phase_Space(REAL sqrt_shat);
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
