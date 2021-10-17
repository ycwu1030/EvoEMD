#ifndef _BASE_MODEL_H_
#define _BASE_MODEL_H_

#include <map>
#include <set>

#include "Neutrino.h"
#include "Particles.h"

class BaseModel {
public:
    typedef std::map<int, Particle *> Particle_List;
    static BaseModel &Get_Particle_Model();
    void Register_Particle(Particle *);
    bool Register_POI(int PID);
    bool Set_Mass(int PID, double mass);

private:
    BaseModel();
    ~BaseModel();

    Particle_List PL;   // All Particle
    std::set<int> POI;  // Particle of Interested
};

#endif  //_BASE_MODEL_H_
