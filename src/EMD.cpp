#include <iostream>
#include "EMD.h"
#include "EffDOF.h"
#include "Physics_Constants.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_sf_bessel.h"

EMD::EMD()
{
    Set_Temperature(1e14,1e5);
}
EMD::EMD(REAL _Ti, REAL _Tr, REAL _Tf, REAL _TInflation)
{
    Set_Temperature(_Ti,_Tr,_Tf,_TInflation);
}
void EMD::Set_Temperature(REAL _Ti, REAL _Tr, REAL _Tf, REAL _TInflation)
{
    Ti = _Ti;
    Tr = _Tr;
    Tf = _Tf;
    TInflation = _TInflation;

    gei = ge(Ti);
    gsi = gs(Ti);
    ger = ge(Tr);
    gsr = gs(Tr);

    Delta = 3.0/4.0*gei/gsi*Ti/Tr;
    Get_Te();

    gee = ge(Te);
    gse = gs(Te);

    CoverD = 5.0/2.0;//BRM;
    Hubble_RD_at_Tr = M_PI*sqrt(ger)/3.0/sqrt(10.0)*Tr*Tr/PHY_MP;

}

double Equation_For_LogTe(double logTe, void *params)
{
    EMD *ml = (EMD*)params;
    REAL Tr = ml->Tr;
    REAL Delta = ml->Delta;
    REAL ger = ml->ger;
    REAL gee = ge(pow(10,logTe));
    REAL gse = gs(pow(10,logTe));
    return 5.0*logTe - 5.0*log10(Tr) - log10(Delta) - log10(4.0/3.0) - log10(gse) - log10(ger) + 2.0*log10(gee);
}

void EMD::Get_Te()
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    // long double logTemax = log10(Tr*pow(Delta*4.0l/3.0l*gsi/ger,1.0l/5.0l))+0.01;
    // long double logTemin = log10(Tr*pow(Delta*4.0l/3.0l*gsr/ger,1.0l/5.0l))-0.01;
    long double logTemax = log10(Tr*pow(Delta*4.0l/3.0l,1.0l/5.0l))+0.01;
    long double logTemin = log10(Tr*pow(Delta*4.0l/3.0l*gsr*ger/gei/gei,1.0l/5.0l))-0.01;
    double r;
    double logTe;
    gsl_function F;
    F.function = &Equation_For_LogTe;
    F.params = this;
    gsl_root_fsolver_set(s, &F, logTemin,logTemax);
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        logTemin = gsl_root_fsolver_x_lower(s);
        logTemax = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(logTemin,logTemax,1e-4,1e-3);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    if (status == GSL_SUCCESS)
    {
        logTe = gsl_root_fsolver_root(s);
    }
    else
    {
        std::cout<<"Not find solution for Te, using dof at Tr to calculate it instead"<<std::endl;
        logTe = log10(Tr*pow(Delta*4.0l/3.0l*gsr/ger,1.0l/5.0l));
    }
    gsl_root_fsolver_free(s);
    
    Te = pow(10,logTe);
}

REAL EMD::Get_Hubble_at_T(REAL Temp)
{
    double geT = ge(Temp);
    double gsT = gs(Temp);

    if (Temp <= Tr)
    {
        return M_PI*sqrt(geT)/3.0/sqrt(10.0)*Temp*Temp/PHY_MP;
    }
    
    if (Temp <= Te)
    {
        return Hubble_RD_at_Tr*geT/ger*pow(Temp/Tr,4.0);
    }

    if (Temp <= Ti)
    {
        return Hubble_RD_at_Tr*sqrt(Delta*4.0/3.0*gsT/ger)*pow(Temp/Tr,1.5);
    }

    return M_PI*sqrt(geT)/3.0/sqrt(10.0)*Temp*Temp/PHY_MP;

}

REAL Number_Density_Eq(REAL T, REAL M, REAL g)
{
    double z = M/T;
    double z2K2 = z*z*gsl_sf_bessel_Kn(2,z);
    return g*pow(T,3)/2/M_PI/M_PI*z2K2;
}

REAL Entropy_Density(REAL T)
{
    return 2*M_PI*M_PI/45.0*106.5*pow(T,3);
}

REAL Yield_Eq(REAL T, REAL M, REAL g)
{
    return Number_Density_Eq(T,M,g)/Entropy_Density(T);
}
