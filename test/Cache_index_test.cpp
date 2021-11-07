#include <bitset>

#include "EvoEMD/Cache.h"

using namespace std;
using namespace EvoEMD;
int main(int argc, char const *argv[]) {
    REAL T = 2.4543343252345345639856739845768934586794875298345982734859234e6;
    INDEX ind_ori = OBTAIN_KEY(T, sizeof(REAL) * 8);
    INDEX ind_masked = OBTAIN_KEY(T);

    cout << ind_ori << " " << bitset<sizeof(INDEX) * 8>(ind_ori) << endl;
    cout << ind_masked << " " << bitset<sizeof(INDEX) * 8>(ind_masked) << endl;

    CACHE cc;
    INDEX idT = OBTAIN_KEY(T);
    cc.Insert(idT, 2.01023);
    REAL T1 = 2.3945829384589279834e-3;
    INDEX idT1 = OBTAIN_KEY(T1);
    cc.Insert(idT1, 3.10495);
    cc.Print_Cache();

    return 0;
}
