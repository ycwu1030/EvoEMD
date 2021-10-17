#ifndef _BASE_MODEL_H_
#define _BASE_MODEL_H_

#include <map>

#include "Neutrino.h"
#include "Particles.h"

class BaseModel {
public:
    typedef std::map<int, Particle *> Particle_List;
    static BaseModel *Get_Model();
    void Register_Particle(Particle *);
    void Register_POI(Particle *);

private:
    BaseModel();
    ~BaseModel();

    Particle_List PL;   // All Particle
    Particle_List POI;  // Particle of Interested

    void Register_Particle(Particle *, Particle_List &pp);
};

#endif  //_BASE_MODEL_H_
