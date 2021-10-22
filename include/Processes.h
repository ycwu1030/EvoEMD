#ifndef _PROCESSES_H_
#define _PROCESSES_H_

#include "CollisionRate.h"
#include "ParticleFactory.h"

template <int ID>
class Amp_N_LPhi : public Amplitude {
public:
    Amp_N_LPhi() {
        N_INITIAL = 1;
        N_FINAL = 2;
        Particle_Factory &pf = Particle_Factory::Get_Particle_Factory();
        INITIAL.push_back(pf.Get_Particle(static_cast<Particle_Factory::Particle_Name>(900000 + ID)));
        FINAL.push_back(pf.Get_Particle(Particle_Factory::Phi));
        FINAL.push_back(pf.Get_Particle(Particle_Factory::DeltaL));
    };
    virtual void Update_Amp(REAL sqrt_shat) { amp_res.n_diag[Process_Amp::SS] = 1; }
    virtual REAL Get_Coeff(REAL T) { return 1; }
};

#endif  //_PROCESSES_H
