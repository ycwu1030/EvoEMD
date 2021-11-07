#ifndef _EVOEMD_CACHE_H_
#define _EVOEMD_CACHE_H_

#include <unordered_map>

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

INDEX OBTAIN_KEY(REAL T, int bit_to_be_compare = 62);

class CACHE {
public:
    typedef std::unordered_map<INDEX, REAL> CACHE_MAP;
    CACHE() = default;
    ~CACHE() = default;
    void Insert(INDEX, REAL);
    bool Get(INDEX, REAL &);
    void Clean_Cache() { cache_data.clear(); }
    void Print_Cache();

private:
    CACHE_MAP cache_data;
};

// inline bool Lookup_Cache(int process_id, REAL T, REAL &res) {
//     CACHE &ca = CACHE::Get_Cache();
//     INDEX id = OBTAIN_KEY(T);
//     return ca.Get(process_id, id, res);
// }

// inline void Make_Cache(int process_id, REAL T, REAL res) {
//     CACHE &ca = CACHE::Get_Cache();
//     INDEX id = OBTAIN_KEY(T);
//     ca.Insert(process_id, id, res);
// }

// inline void Clean_Cache() { CACHE::Get_Cache().Clean_Cache(); }

}  // namespace EvoEMD

#endif  //_EVOEMD_CACHE_H_
