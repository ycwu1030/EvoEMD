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

Pseudo_Particle *Particle_Factory::Get_Particle(Particle_Name pn) {
    Particle_List::iterator iter = PL.find(pn);
    if (iter == PL.end()) {
        Pseudo_Particle *tmp = Register_Default_Particle(pn);
        return tmp;
    } else {
        return iter->second;
    }
}

Pseudo_Particle *Particle_Factory::Register_Default_Particle(Particle_Name pn) {
    Pseudo_Particle *tmp = nullptr;
    switch (pn) {
        case RHN1:
            tmp = new Fermion(1000.0, static_cast<int>(pn), 2, true);
            break;
        case RHN2:
            tmp = new Fermion(2000.0, static_cast<int>(pn), 2, true);
            break;
        case RHN3:
            tmp = new Fermion(3000.0, static_cast<int>(pn), 2, true);
            break;
        case DeltaL:
            tmp = new Fermion(static_cast<int>(pn), 3 * (2 + 1));
            break;
        case LNU1:
            tmp = new Fermion(0.05 * eV, static_cast<int>(pn), 2);
            break;
        case Phi:
            tmp = new Boson(static_cast<int>(pn), 2);
            break;
        case Chi:
            tmp = new Fermion(100.0, static_cast<int>(pn), 2);
            break;
        case S:
            tmp = new Boson(200.0, static_cast<int>(pn), 1);
            break;
        default:
            break;
    }
    if (!tmp) {
        PL.insert(std::make_pair(pn, tmp));
    }
    return tmp;
}

bool Particle_Factory::Register_POI(Particle_Name pn) {
    Particle_List::iterator iter = PL.find(pn);
    if (iter == PL.end()) {
        Pseudo_Particle *tmp = Register_Default_Particle(pn);
        if (!tmp) {
            POI.insert(pn);
            return true;
        } else {
            return false;
        }
    } else {
        POI.insert(pn);
        return true;
    }
}

bool Particle_Factory::Set_Mass(Particle_Name PID, double mass) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        return false;
    } else {
        iter->second->Set_Mass(mass);
        return true;
    }
}
