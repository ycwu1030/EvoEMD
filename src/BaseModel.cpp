#include "BaseModel.h"

BaseModel::BaseModel() {}

BaseModel::~BaseModel() {
    Particle_List::iterator iter;
    for (iter = PL.begin(); iter != PL.end(); ++iter) {
        delete iter->second;
    }
}

BaseModel &BaseModel::Get_Particle_Model() {
    static BaseModel BM;
    return BM;
}

void BaseModel::Register_Particle(Particle *part) {
    int PID_key = part->Get_PID();
    Particle_List::iterator iter = PL.find(PID_key);
    if (iter == PL.end()) {
        PL.insert(std::make_pair(PID_key, part));
    } else {
        iter->second = part;
    }
}

bool BaseModel::Register_POI(int PID) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        return false;
    } else {
        POI.insert(PID);
    }
}

bool BaseModel::Set_Mass(int PID, double mass) {
    Particle_List::iterator iter = PL.find(PID);
    if (iter == PL.end()) {
        return false;
    } else {
        iter->second->Set_Mass(mass);
        return true;
    }
}
