#ifndef _EVOEMD_CACHE_H_
#define _EVOEMD_CACHE_H_

#include <unordered_map>

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

INDEX OBTAIN_KEY(REAL T, int bit_to_be_compare = 62);

class CACHE {
public:
    typedef std::unordered_map<INDEX, REAL> CACHE_MAP;
    static CACHE &Get_Cache();
    void Insert(INDEX, REAL);
    bool Get(INDEX, REAL &);
    void Clean_Cache() { cache_data.clear(); }
    void Print_Cache();

private:
    CACHE() = default;
    ~CACHE() = default;
    CACHE(const CACHE &) = delete;
    CACHE(CACHE &&) = delete;
    CACHE_MAP cache_data;

    CACHE &operator=(const CACHE &) = delete;
    CACHE &operator=(CACHE &&) = delete;
};

inline bool Lookup_Cache(REAL T, REAL &res) {
    CACHE &ca = CACHE::Get_Cache();
    INDEX id = OBTAIN_KEY(T);
    return ca.Get(id, res);
}

inline void Make_Cache(REAL T, REAL res) {
    CACHE &ca = CACHE::Get_Cache();
    INDEX id = OBTAIN_KEY(T);
    ca.Insert(id, res);
}

inline void Clean_Cache() { CACHE::Get_Cache().Clean_Cache(); }

}  // namespace EvoEMD

#endif  //_EVOEMD_CACHE_H_
