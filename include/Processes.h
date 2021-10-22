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
    virtual void Update_Amp(REAL sqrt_shat) {
        amp_res.n_diag = 1;

        using PS = Process_Amp::PROPAGATOR_STRUCTURE;
        using CRF = Process_Amp::CTH_RES_FULL;
        // * For decay, no propagator, using dummy propagator;
        CRF dummy_CRF = {1};

        PS dummy_propagator = std::make_pair(0, dummy_CRF);
        amp_res.amps_denominator.clear();
        amp_res.amps_denominator.push_back(std::make_pair(dummy_propagator, dummy_propagator));

        // * Numerator:
        CRF num(3);
        Nu_TypeI_SeeSaw &nm = Nu_TypeI_SeeSaw::Get_Neutrino_Model();
        num[0] = 2 * 2 * (INITIAL[0]->Get_Mass()) * nm.Get_YdagYij(ID, ID);
        num[1] = 0;
        num[2] = 0;
    }
    virtual REAL Get_Coeff(REAL T) { return 1; }
};

#endif  //_PROCESSES_H
