#ifndef _PARTICLE_FACTORY_H_
#define _PARTICLE_FACTORY_H_

#include <map>
#include <set>

#include "Neutrino.h"
#include "Particles.h"

class Particle_Factory {
public:
    typedef enum {
        RHN1 = 900001,
        RHN2 = 900002,
        RHN3 = 900003,
        DeltaL = 900011,
        Phi = 900025,
        Chi = 900100,
        S = 900200
    } Particle_Name;
    typedef std::map<Particle_Name, Pseudo_Particle *> Particle_List;
    static Particle_Factory &Get_Particle_Factory();
    Pseudo_Particle *Get_Particle(Particle_Name);
    bool Register_POI(Particle_Name);
    bool Set_Mass(Particle_Name PID, double mass);

private:
    Particle_Factory();
    ~Particle_Factory();

    Particle_List PL;             // All Particle
    std::set<Particle_Name> POI;  // Particle of Interested

    Pseudo_Particle *Register_Default_Particle(Particle_Name);
};

#endif  //_PARTICLE_FACTORY_H_
