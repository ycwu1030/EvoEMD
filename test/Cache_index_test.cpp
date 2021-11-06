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

    Make_Cache(T, 1.101);
    Make_Cache(T, 2.01023);
    REAL T1 = 2.3945829384589279834e-3;
    Make_Cache(T1, 3.10495);
    CACHE::Get_Cache().Print_Cache();

    return 0;
}
