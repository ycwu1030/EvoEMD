#include "EMD.h"

#include <iostream>

#include "EffDOF.h"
#include "Physics_Constants.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_zeta.h"
#include "spdlog_wrapper.h"

EMD::EMD() {
    SPDLOG_INFO_FILE("DEFAULT CONSTRUCTOR FOR EMD.");
    Set_Temperature(1e14, 1e5);
}
EMD::EMD(REAL _Ti, REAL _Tr, REAL _Tf, REAL _TInflation) {
    SPDLOG_INFO_FILE("CONSTRUCTOR FOR EMD with Ti = {:+9.8e}, Tr = {:+9.8e}, Tf = {:+9.8e}, Tinflation = {:+9.8e}", _Ti,
                     _Tr, _Tf, _TInflation);
    Set_Temperature(_Ti, _Tr, _Tf, _TInflation);
}
void EMD::Set_Temperature(REAL _Ti, REAL _Tr, REAL _Tf, REAL _TInflation) {
    SPDLOG_INFO_FILE(
        "Setting Temperatures for EMD with Ti = {:+9.8e}, Tr = {:+9.8e}, Tf = {:+9.8e}, Tinflation = {:+9.8e}", _Ti,
        _Tr, _Tf, _TInflation);
    Ti = _Ti;
    Tr = _Tr;
    Tf = _Tf;
    TInflation = _TInflation;

    gei = ge(Ti);
    gsi = gs(Ti);
    ger = ge(Tr);
    gsr = gs(Tr);

    Delta = 3.0 / 4.0 * gei / gsi * Ti / Tr;
    Solve_Te();

    SPDLOG_INFO_FILE("Setting Temperature Te = {:+9.8e}.", Te);
    gee = ge(Te);
    gse = gs(Te);

    CoverD = 5.0 / 2.0;  // BRM;
    Hubble_RD_at_Tr = M_PI * sqrt(ger) / 3.0 / sqrt(10.0) * Tr * Tr / PHY_MP;
}

double Equation_For_LogTe(double logTe, void *params) {
    EMD *ml = (EMD *)params;
    REAL Tr = ml->Tr;
    REAL Delta = ml->Delta;
    REAL ger = ml->ger;
    REAL gee = ge(pow(10, logTe));
    REAL gse = gs(pow(10, logTe));
    return 5.0 * logTe - 5.0 * log10(Tr) - log10(Delta) - log10(4.0 / 3.0) - log10(gse) - log10(ger) + 2.0 * log10(gee);
}

void EMD::Solve_Te() {
    SPDLOG_DEBUG_FILE("Try to Get Te from Ti ({:+9.8e}) and Tr ({:+9.8e})", Ti, Tr);
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    // long double logTemax =
    // log10(Tr*pow(Delta*4.0l/3.0l*gsi/ger,1.0l/5.0l))+0.01; long double logTemin
    // = log10(Tr*pow(Delta*4.0l/3.0l*gsr/ger,1.0l/5.0l))-0.01;
    long double logTemax = log10(Tr * pow(Delta * 4.0l / 3.0l, 1.0l / 5.0l)) + 0.01;
    long double logTemin = log10(Tr * pow(Delta * 4.0l / 3.0l * gsr * ger / gei / gei, 1.0l / 5.0l)) - 0.01;
    SPDLOG_DEBUG_FILE("The range for finding Te is [{:+9.8e},{:+9.8e}]", pow(10, logTemin), pow(10, logTemax));
    double r;
    double logTe;
    gsl_function F;
    F.function = &Equation_For_LogTe;
    F.params = this;
    gsl_root_fsolver_set(s, &F, logTemin, logTemax);
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        logTemin = gsl_root_fsolver_x_lower(s);
        logTemax = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(logTemin, logTemax, 1e-4, 1e-3);
        SPDLOG_DEBUG_FILE("At Iter = {}, the interval is [{:+9.8e},{:+9.8e}]", iter, pow(10, logTemin),
                          pow(10, logTemax));
    } while (status == GSL_CONTINUE && iter < max_iter);

    if (status == GSL_SUCCESS) {
        logTe = gsl_root_fsolver_root(s);
    } else {
        logTe = log10(Tr * pow(Delta * 4.0l / 3.0l * gsr / ger, 1.0l / 5.0l));
        SPDLOG_WARN_FILE("NOT FINDING SOLUTION FOR Te, USING NAIVE ESTIMATIONS: Te = {:+9.8e}.", pow(10, logTe));
    }
    gsl_root_fsolver_free(s);

    Te = pow(10, logTe);
    SPDLOG_DEBUG_FILE("Solution is Te = {:+9.8e}.", Te);
}

REAL EMD::Get_Hubble_at_T(REAL Temp) {
    SPDLOG_INFO_FILE("Calculate Hubble at T = {:+9.8e}.", Temp);
    REAL geT = ge(Temp);
    REAL gsT = gs(Temp);
    REAL HatT;
    if (Temp <= Tr) {
        HatT = M_PI * sqrt(geT) / 3.0 / sqrt(10.0) * Temp * Temp / PHY_MP;
        SPDLOG_INFO_FILE("Late Radiation Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
        return HatT;
    }

    if (Temp <= Te) {
        HatT = Hubble_RD_at_Tr * geT / ger * pow(Temp / Tr, 4.0);
        SPDLOG_INFO_FILE("Entropy Production Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
        return HatT;
    }

    if (Temp <= Ti) {
        HatT = Hubble_RD_at_Tr * sqrt(Delta * 4.0 / 3.0 * gsT / ger) * pow(Temp / Tr, 1.5);
        SPDLOG_INFO_FILE("Early Matter Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
        return HatT;
    }

    HatT = M_PI * sqrt(geT) / 3.0 / sqrt(10.0) * Temp * Temp / PHY_MP;
    SPDLOG_INFO_FILE("Early Radiation Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
    return HatT;
}

REAL EMD::Get_Hubble_at_T_Period(Period pd, REAL Temp) {
    SPDLOG_INFO_FILE("Calculate Hubble For Period {} at T = {:+9.8e}", pd, Temp);
    REAL geT = ge(Temp);
    REAL gsT = gs(Temp);
    REAL HatT;
    switch (pd) {
        case ERDE:
            HatT = M_PI * sqrt(geT) / 3.0 / sqrt(10.0) * Temp * Temp / PHY_MP;
            SPDLOG_INFO_FILE("Early Radiation Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
            return HatT;
        case EMDE:
            HatT = Hubble_RD_at_Tr * sqrt(Delta * 4.0 / 3.0 * gsT / ger) * pow(Temp / Tr, 1.5);
            SPDLOG_INFO_FILE("Early Matter Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
            return HatT;
        case EPE:
            HatT = Hubble_RD_at_Tr * geT / ger * pow(Temp / Tr, 4.0);
            SPDLOG_INFO_FILE("Entropy Production Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
            return HatT;
        case RDE:
            HatT = M_PI * sqrt(geT) / 3.0 / sqrt(10.0) * Temp * Temp / PHY_MP;
            SPDLOG_INFO_FILE("Late Radiation Era, H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
            return HatT;
        default:
            HatT = M_PI * sqrt(geT) / 3.0 / sqrt(10.0) * Temp * Temp / PHY_MP;
            SPDLOG_WARN_FILE("No Matching Period, Using RD instead: H = {:+9.8e}, T = {:+9.8e}.", HatT, Temp);
            return HatT;
    }
}

REAL Number_Density_Eq(REAL T, REAL M, REAL g) {
    REAL z = M / T;
    SPDLOG_INFO_FILE(
        "Calculating Particle Number Density at Equilibrium, T = {:+9.8e}, M = {:+9.8e}, g = {}, z = {:+9.8e}", T, M, g,
        z);
    REAL z2K2 = z * z * gsl_sf_bessel_Kn(2, z);
    SPDLOG_DEBUG_FILE("For Calculating Particle Number Density at Eq, z^2*K2@z = {:+9.8e}", z2K2);
    REAL neq = g * pow(T, 3) / 2 / M_PI / M_PI * z2K2;
    SPDLOG_DEBUG_FILE("For Calculating Particle Number Density at Eq, neq = {:+9.8e}", neq);
    return neq;
}

REAL Number_Density_Eq_BE(REAL T, REAL g) {
    SPDLOG_INFO_FILE("Calculating Particle Number Density at Equilibrium for Massless Boson, T = {:+9.8e}, g = {}", T,
                     g);
    REAL neq = g * pow(T, 3) / M_PI / M_PI;
    SPDLOG_DEBUG_FILE("For Calculating Particle Number Density at Eq, neq = {:+9.8e}", neq);
    return neq;
}

REAL Number_Density_Eq_FD(REAL T, REAL g) {
    SPDLOG_INFO_FILE("Calculating Particle Number Density at Equilibrium for Massless Fermion, T = {:+9.8e}, g = {}", T,
                     g);
    REAL neq = 3.0 / 4.0 * g * gsl_sf_zeta_int(3.0) * pow(T, 3) / M_PI / M_PI;
    SPDLOG_DEBUG_FILE("For Calculating Particle Number Density at Eq, neq = {:+9.8e}", neq);
    return neq;
}

REAL Entropy_Density(REAL T) {
    REAL gsT = gs(T);

    REAL s = 2 * M_PI * M_PI / 45.0 * gsT * pow(T, 3);
    SPDLOG_DEBUG_FILE("Entropy density s = {:+9.8e} at T = {:+9.8e}.", s, T);
    return s;
}

REAL Yield_Eq(REAL T, REAL M, REAL g) {
    REAL yeq = Number_Density_Eq(T, M, g) / Entropy_Density(T);
    SPDLOG_DEBUG_FILE("Yield at Equilibrium Yeq = {:+9.8e} at T = {:+9.8e}, for M = {:+9.8e}, g = {}.", yeq, T, M, g);
    return yeq;
}

REAL Yield_Eq_BE(REAL T, REAL g) {
    REAL yeq = Number_Density_Eq_BE(T, g) / Entropy_Density(T);
    SPDLOG_DEBUG_FILE("Yield at Equilibrium for Massless Boson Yeq = {:+9.8e} at T = {:+9.8e}, g = {}.", yeq, T, g);
    return yeq;
}

REAL Yield_Eq_FD(REAL T, REAL g) {
    REAL yeq = Number_Density_Eq_FD(T, g) / Entropy_Density(T);
    SPDLOG_DEBUG_FILE("Yield at Equilibrium for Massless Fermion Yeq = {:+9.8e} at T = {:+9.8e}, g = {}.", yeq, T, g);
    return yeq;
}
