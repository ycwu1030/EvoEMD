#include <cmath>
#include <fstream>
#include <iostream>

#include "EvoEMD/Hubble.h"

using namespace std;
using namespace EvoEMD;
int main(int argc, char const *argv[]) {
    Hubble_History EMDE(1e14, 10);
    ofstream output("Hubble_History.txt");
    output << "T\tH" << endl;
    for (REAL lT = 14.5; lT >= 0; lT -= 0.1) {
        REAL T = pow(10, lT);
        output << T << "\t" << EMDE.Get_Hubble_at_T(T) << endl;
    }

    return 0;
}
