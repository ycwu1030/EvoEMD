#ifndef _PARTICLE_FACTORY_H_
#define _PARTICLE_FACTORY_H_

#include <map>
#include <set>

#include "EvoEMD/Neutrino.h"
#include "EvoEMD/ParticlesBase.h"

namespace EvoEMD {
class Particle_Factory {
public:
    typedef std::map<int, Pseudo_Particle *> Particle_List;

    static Particle_Factory &Get_Particle_Factory();

    Pseudo_Particle *Get_Particle(int PID);
    bool Register_Particle(Pseudo_Particle *);
    bool Register_POI(int PID);
    bool Set_Mass(int PID, double mass);

private:
    Particle_Factory();
    ~Particle_Factory();

    Particle_List PL;   // All Particle
    std::set<int> POI;  // Particle of Interested
};
}  // namespace EvoEMD
#endif  //_PARTICLE_FACTORY_H_
