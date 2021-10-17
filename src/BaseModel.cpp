#include "BaseModel.h"

BaseModel::BaseModel() {}

BaseModel::~BaseModel() {
    Particle_List::iterator iter;
    for (iter = PL.begin(); iter != PL.end(); ++iter) {
        delete iter->second;
    }
}

BaseModel *BaseModel::Get_Model() {
    static BaseModel BM;
    return &BM;
}

void BaseModel::Register_Particle(Particle *part, Particle_List &pp) {
    int PID_key = part->Get_PID();
    Particle_List::iterator iter = pp.find(PID_key);
    if (iter == pp.end()) {
        pp.insert(std::make_pair(PID_key, part));
    } else {
        iter->second = part;
    }
}

void BaseModel::Register_Particle(Particle *part) { Register_Particle(part, PL); }
void BaseModel::Register_POI(Particle *part) { Register_Particle(part, POI); }
