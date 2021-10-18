#ifndef _COLLISION_RATE_H_
#define _COLLISION_RATE_H_

#include "Particles.h"
#include "RealTypes.h"

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

public:
    Collision_Rate(){};
    virtual ~Collision_Rate(){};

    virtual REAL Get_CP_Conserving_Rate(REAL T) = 0;
    virtual REAL Get_CP_Violating_Rate(REAL T) = 0;
};

class Process {
public:
    typedef std::vector<Pseudo_Particle *> INITIAL_PARTICLES;
    typedef std::vector<Pseudo_Particle *> FINAL_PARTICLES;

    Process();
    ~Process();

protected:
    INITIAL_PARTICLES INIT;
    FINAL_PARTICLES FINAL;

    virtual REAL Matrix_Element_Square_Integrated_over_Final_State_Phase_Space(REAL T) = 0;
};

class Decay12_Rate : public Collision_Rate {
    // * For two body decay
private:
public:
    Decay12_Rate();
    ~Decay12_Rate();

    virtual REAL Get_CP_Conserving_Rate(REAL T);
    virtual REAL Get_CP_Violating_Rate(REAL T);
};

class Scatter22_Rate {
protected:
public:
    Scatter22_Rate();
    ~Scatter22_Rate();

    virtual REAL Get_CP_Conserving_Rate(REAL T);
    virtual REAL Get_CP_Violating_Rate(REAL T);
};

#endif  //_COLLISION_RATE_H_
