#ifndef _HUBBLE_H_
#define _HUBBLE_H_

#include "EvoEMD/ParameterBase.h"
#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

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
    Hubble_For_Single_Period(const REAL T_start, const REAL T_end, const bool Isentropic = true,
                             const double beta_T = 1, const double beta_s = 3);
    virtual ~Hubble_For_Single_Period(){};

    /**
     * @brief Calculate the Hubble parameter at temperature T, need to be overwrite by derived class
     * @note
     * @param  T: the temperature
     * @retval The hubble parameter
     */
    virtual REAL Get_Hubble_at_T(const REAL T) = 0;
    static REAL Get_Hubble_For_RD(const REAL T);
    REAL Get_T_Start() const { return T_start; }
    REAL Get_T_End() const { return T_end; }
    double Get_beta_T() const { return beta_T; }
    double Get_beta_s() const { return beta_s; }
    bool Is_Isentropic() const { return Isentropic; }
    void Print() const;
};

class Hubble_RD : public Hubble_For_Single_Period {
public:
    Hubble_RD(const REAL T_start, const REAL T_end);
    ~Hubble_RD(){};

    virtual REAL Get_Hubble_at_T(const REAL T) override;
};

class Hubble_EMD : public Hubble_For_Single_Period {
private:
    REAL HRD_at_T_start;

public:
    Hubble_EMD(const REAL T_start, const REAL T_end);
    ~Hubble_EMD(){};

    virtual REAL Get_Hubble_at_T(const REAL T) override;
};

class Hubble_EP : public Hubble_For_Single_Period {
private:
    REAL HRD_at_T_end;
    REAL ge_at_T_end;

public:
    Hubble_EP(const REAL T_start, const REAL T_end);
    ~Hubble_EP(){};

    virtual REAL Get_Hubble_at_T(const REAL T) override;
};

class Hubble_History : public Parameter_Base {
private:
    std::vector<Hubble_For_Single_Period *> Periods;
    std::vector<REAL> Temperatures;
    REAL TRH;
    REAL Ti;
    REAL Te;
    REAL Tr;
    REAL Tf;

    void Solve_Te();

    Hubble_History();
    Hubble_History(const Hubble_History &HH);
    Hubble_History &operator=(const Hubble_History &HH);
    ~Hubble_History();

public:
    static Hubble_History &Get_Hubble_History();
    int Get_N_Period() const { return Periods.size(); }
    int Get_Period_ID_at_T(const REAL T);
    REAL Get_Hubble_at_T(const REAL T);
    double Get_beta_T_at_T(const REAL T);
    Hubble_For_Single_Period *operator[](const int pid);
    Hubble_For_Single_Period *at(const int pid);

    virtual void Update_Value(REAL input) override;
    void Print_History();
};

}  // namespace EvoEMD

#endif  //_HUBBLE_H_
