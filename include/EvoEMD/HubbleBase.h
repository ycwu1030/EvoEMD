#ifndef _HUBBLE_H_
#define _HUBBLE_H_

#include "EvoEMD/HubbleEvolution.h"
#include "EvoEMD/ParameterBase.h"
#include "EvoEMD/RealTypes.h"
#include "gsl/gsl_spline.h"

namespace EvoEMD {

REAL Get_Hubble_For_RD_at_T(REAL T);

class Hubble_Base {
public:
    Hubble_Base(){};
    virtual ~Hubble_Base(){};

    virtual REAL Get_Hubble_at_T(const REAL T) = 0;
    virtual REAL Get_dlna_dlnT_at_T(const REAL T) = 0;
    virtual void Print(){};
};

class Hubble_Splitting_Period : public Hubble_Base {
protected:
    REAL T_BEGIN;
    REAL T_END;
    double beta_R;

public:
    Hubble_Splitting_Period(const REAL T_begin, const REAL T_end, const double beta_R = 1);
    virtual ~Hubble_Splitting_Period(){};

    REAL Get_T_Begin() const { return T_BEGIN; }
    REAL Get_T_End() const { return T_END; }
    double Get_beta_R() const { return beta_R; }

    virtual REAL Get_Hubble_at_T(const REAL T) override;
    virtual REAL Get_dlna_dlnT_at_T(const REAL T) override final;
    virtual void Print() override;
};

class Hubble_RD : public Hubble_Splitting_Period {
public:
    Hubble_RD(const REAL T_start, const REAL T_end);
    ~Hubble_RD(){};
};

class Hubble_EMD : public Hubble_Splitting_Period {
private:
    REAL HRD_at_T_Begin;
    REAL gs_at_T_Begin;

public:
    Hubble_EMD(const REAL T_start, const REAL T_end);
    ~Hubble_EMD(){};

    virtual REAL Get_Hubble_at_T(const REAL T) override final;
};

class Hubble_EP : public Hubble_Splitting_Period {
private:
    REAL HRD_at_T_end;
    REAL ge_at_T_end;

public:
    Hubble_EP(const REAL T_start, const REAL T_end);
    ~Hubble_EP(){};

    virtual REAL Get_Hubble_at_T(const REAL T) override final;
};

class Hubble_Splitting : public Parameter_Base, public Hubble_Base {
private:
    std::vector<Hubble_Base *> Periods;
    std::vector<REAL> Temperatures;
    REAL TRH;
    REAL Ti;
    REAL Te;
    REAL Tr;
    REAL Tf;

    void Solve_Te();

public:
    Hubble_Splitting();
    Hubble_Splitting(const Hubble_Splitting &HH);
    Hubble_Splitting &operator=(const Hubble_Splitting &HH);
    ~Hubble_Splitting();
    // static Hubble_Base *Get_Hubble_Calculator();
    int Get_N_Period() const { return Periods.size(); }
    int Get_Period_ID_at_T(const REAL T);
    Hubble_Base *operator[](const int pid);
    Hubble_Base *at(const int pid);

    virtual REAL Get_Hubble_at_T(const REAL T) override;
    virtual REAL Get_dlna_dlnT_at_T(const REAL T) override;
    virtual void Update_Value(REAL input) override;
    virtual void Print() override;
};

class Hubble_BE : public Parameter_Base, public Hubble_Base {
private:
    VD List_T;
    VD List_U;  // U = k*a
    VD List_H;
    VD List_lnT;
    VD List_lnU;
    VD List_lnH;
    VD List_rhoR;
    VD List_rhoM;

    gsl_interp_accel *acc_Us;
    gsl_interp_accel *acc_Hs;
    gsl_spline *spline_Us;
    gsl_spline *spline_Hs;

    Hubble_Evolution HE;
    REAL Calc_Tr(VD U, VD Y1, VD Y2, REAL krhoR);
    void Clean();

public:
    Hubble_BE();
    ~Hubble_BE(){};
    // static Hubble_Base *Get_Hubble_Calculator();
    virtual REAL Get_Hubble_at_T(const REAL T) override;
    virtual REAL Get_dlna_dlnT_at_T(const REAL T) override;
    virtual void Print() override;
    virtual void Update_Value(REAL input) override;
    VD Get_T_List();
    VD Get_U_List();
    VD Get_H_List();
    VD Get_rhoR_List();
    VD Get_rhoM_List();
};

class Hubble_Factory {
public:
    typedef std::map<unsigned int, Hubble_Base *> Hubble_Calculator_List;

    static Hubble_Factory &Get_Hubble_Factory();
    static Hubble_Base *Get_Hubble_Calculator(int id = -1);
    static void Register_Hubble_Calculator(int id, Hubble_Base *);

    Hubble_Base *Get_Calculator(int id = -1);
    void Register_Calculator(int id, Hubble_Base *);

private:
    Hubble_Factory();
    ~Hubble_Factory();

    Hubble_Calculator_List HCL;
};

}  // namespace EvoEMD

#endif  //_HUBBLE_H_
