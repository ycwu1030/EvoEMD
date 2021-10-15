#ifndef _HUBBLE_H_
#define _HUBBLE_H_

#include "RealTypes.h"

class Hubble_For_Single_Period {
protected:
    REAL T_start;
    REAL T_end;
    bool Isentropic;
    double beta_T;
    double beta_s;

public:
    /**
     * @brief The contructor, default is for a RD or EMD period
     * @note
     * @param  T_start: the starting temperature of this period
     * @param  T_end: the ending temperature of this period
     * @param  Isentropic: whether this period is a isentropic period.
     * @param  beta_T: T~a^{-beta_T}
     * @param  beta_s: s~a^{-beta_s} S=s a^3
     * @retval
     */
    Hubble_For_Single_Period(REAL T_start, REAL T_end, bool Isentropic = true, double beta_T = 1, double beta_s = 3);
    virtual ~Hubble_For_Single_Period(){};

    /**
     * @brief Calculate the Hubble parameter at temperature T, need to be overwrite by derived class
     * @note
     * @param  T: the temperature
     * @retval The hubble parameter
     */
    virtual REAL Get_Hubble_at_T(REAL T) = 0;
    static REAL Get_Hubble_For_RD(REAL T);
    REAL Get_T_Start() const { return T_start; }
    REAL Get_T_End() const { return T_end; }
    double Get_beta_T() const { return beta_T; }
    double Get_beta_s() const { return beta_s; }
    bool Is_Isentropic() const { return Isentropic; }
};

class Hubble_RD : public Hubble_For_Single_Period {
public:
    Hubble_RD(REAL T_start, REAL T_end);
    ~Hubble_RD(){};

    virtual REAL Get_Hubble_at_T(REAL T);
};

class Hubble_EMD : public Hubble_For_Single_Period {
private:
    REAL HRD_at_T_start;

public:
    Hubble_EMD(REAL T_start, REAL T_end);
    ~Hubble_EMD(){};

    virtual REAL Get_Hubble_at_T(REAL T);
};

#endif  //_HUBBLE_H_
