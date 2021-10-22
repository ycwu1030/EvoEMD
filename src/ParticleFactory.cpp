#include "ParticleFactory.h"

Particle_Factory::Particle_Factory() {}

Particle_Factory::~Particle_Factory() {
    Particle_List::iterator iter;
    for (iter = PL.begin(); iter != PL.end(); ++iter) {
        delete iter->second;
    }
}

Particle_Factory &Particle_Factory::Get_Particle_Factory() {
    static Particle_Factory PF;
    return PF;
}

Pseudo_Particle *Particle_Factory::Get_Particle(int PID) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        std::cout << "Cannot find particle with PID = " << PID << std::endl;
        return nullptr;
    } else {
        return iter->second;
    }
}

bool Particle_Factory::Register_POI(int PID) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        std::cout << "Cannot find particle with PID = " << PID << std::endl;
        return false;
    } else {
        POI.insert(PID);
        return true;
    }
}

bool Particle_Factory::Set_Mass(int PID, double mass) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        std::cout << "Cannot find particle with PID = " << PID << std::endl;
        return false;
    } else {
        iter->second->Set_Mass(mass);
        return true;
    }
}
