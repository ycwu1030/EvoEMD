#include "EvoEMD/Cache.h"

namespace EvoEMD {

INDEX OBTAIN_KEY(REAL T, int bit_to_be_compare) {
    int bit_difference = sizeof(INDEX) * 8 - bit_to_be_compare;
    if (bit_difference < 0) {
        bit_difference = 0;
    }
    INDEX ONE = 1;
    // * mask will be:
    // * 11....110..0
    // * with bit_to_be_compare 1's and (sizeof(INDEX)*8-bit_to_be_compare) 0's
    INDEX mask = -(ONE << bit_difference);
    INDEX res = mask & (*(INDEX *)&T);
    return res;
}

void CACHE::Insert(INDEX ind, REAL res) { cache_data[ind] = res; }

bool CACHE::Get(INDEX ind, REAL &res) {
    int exist = cache_data.count(ind);
    if (exist == 0) {
        res = 0;
        return false;
    }
    res = cache_data[ind];
    return true;
}

void CACHE::Print_Cache() {
    for (auto &&x : cache_data) {
        std::cout << x.first << " = " << x.second << std::endl;
    }
}

}  // namespace EvoEMD
