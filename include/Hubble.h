#ifndef _HUBBLE_H_
#define _HUBBLE_H_

#include "RealTypes.h"

class Hubble_For_Single_Period {
private:
    bool Isentropic;
    double beta_T;
    double beta_s;
    REAL Get_Hubble_For_RD(REAL T);

public:
    /**
     * @brief The contructor, default is for a RD or EMD period
     * @note
     * @param  Isentropic: whether this period is a isentropic period.
     * @param  beta_T: T~a^{-beta_T}
     * @param  beta_s: s~a^{-beta_s} S=s a^3
     * @retval
     */
    Hubble_For_Single_Period(bool Isentropic = true, double beta_T = 1, double beta_s = 3);
    virtual ~Hubble_For_Single_Period(){};

    virtual REAL Get_Hubble_at_T(REAL T) = 0;
    double Get_beta_T() const { return beta_T; }
    double Get_beta_s() const { return beta_s; }
    bool Is_Isentropic() const { return Isentropic; }
};

#endif  //_HUBBLE_H_
