#include <iostream>

#include "Particles.h"

using namespace std;

int main(int argc, char const *argv[]) {
    Pseudo_Particle *const ptcl = new Fermion(100.1, 900001, 2);
    typedef set<Pseudo_Particle *const> PS;
    PS ps;
    ps.insert(ptcl);
    cout << ptcl->Get_Mass() << endl;
    PS::iterator iter = ps.begin();
    (*iter)->Set_Mass(1);
    // cout << ptcl->Get_Mass() << endl;
    return 0;
}
