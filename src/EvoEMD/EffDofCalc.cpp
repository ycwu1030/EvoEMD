#include <gsl/gsl_spline.h>

#include <cmath>

#include "EvoEMD/EffDOF.h"
#include "EvoEMD/EffDofTable.inc"

namespace EvoEMD {
class EffDofCalc {
private:
    EffDofCalc();
    ~EffDofCalc();
    REAL T_Points[N_T_gs_ge_Points];
    REAL gs_Points[N_T_gs_ge_Points];
    REAL ge_Points[N_T_gs_ge_Points];
    gsl_interp_accel *gs_acc;
    gsl_interp_accel *ge_acc;
    gsl_spline *gs_spline;
    gsl_spline *ge_spline;
    REAL T_MIN;
    REAL T_MAX;

public:
    static EffDofCalc &Get_Calculator() {
        static EffDofCalc calc;
        return calc;
    };
    REAL Calc_ge(REAL T);
    REAL Calc_gs(REAL T);
    REAL Calc_dlnge_dlnT(REAL T);
    REAL Calc_dlngs_dlnT(REAL T);
};

EffDofCalc::EffDofCalc() {
    for (int i = 0; i < N_T_gs_ge_Points; i++) {
        T_Points[i] = T_gs_ge_Points[i][0];
        gs_Points[i] = T_gs_ge_Points[i][1];
        ge_Points[i] = T_gs_ge_Points[i][2];
    }
    T_MIN = T_Points[0];
    T_MAX = T_Points[N_T_gs_ge_Points - 1];
    gs_acc = gsl_interp_accel_alloc();
    ge_acc = gsl_interp_accel_alloc();
    gs_spline = gsl_spline_alloc(gsl_interp_steffen, N_T_gs_ge_Points);
    ge_spline = gsl_spline_alloc(gsl_interp_steffen, N_T_gs_ge_Points);

    gsl_spline_init(gs_spline, T_Points, gs_Points, N_T_gs_ge_Points);
    gsl_spline_init(ge_spline, T_Points, ge_Points, N_T_gs_ge_Points);
}

EffDofCalc::~EffDofCalc() {
    gsl_spline_free(ge_spline);
    gsl_spline_free(gs_spline);
    gsl_interp_accel_free(ge_acc);
    gsl_interp_accel_free(gs_acc);
}

REAL EffDofCalc::Calc_ge(REAL T) {
    if (T <= T_MIN) return ge_Points[0];
    if (T >= T_MAX) return ge_Points[N_T_gs_ge_Points - 1];
    return gsl_spline_eval(ge_spline, T, ge_acc);
}

REAL EffDofCalc::Calc_gs(REAL T) {
    if (T <= T_MIN) return gs_Points[0];
    if (T >= T_MAX) return gs_Points[N_T_gs_ge_Points - 1];
    return gsl_spline_eval(gs_spline, T, gs_acc);
}

REAL EffDofCalc::Calc_dlnge_dlnT(REAL T) {
    if (T <= T_MIN) return 0;
    if (T >= T_MAX) return 0;
    REAL gge = Calc_ge(T);
    REAL dge_dT = gsl_spline_eval_deriv(ge_spline, T, ge_acc);
    return dge_dT * T / gge;
}

REAL EffDofCalc::Calc_dlngs_dlnT(REAL T) {
    if (T <= T_MIN) return 0;
    if (T >= T_MAX) return 0;
    REAL ggs = Calc_gs(T);
    REAL dgs_dT = gsl_spline_eval_deriv(gs_spline, T, gs_acc);
    return dgs_dT * T / ggs;
}

REAL f_ge(REAL T) { return EffDofCalc::Get_Calculator().Calc_ge(T); }
REAL f_gs(REAL T) { return EffDofCalc::Get_Calculator().Calc_gs(T); }
REAL f_ge_star(REAL T) { return 1.0 + EffDofCalc::Get_Calculator().Calc_dlnge_dlnT(T) / 4.0; }
REAL f_gs_star(REAL T) { return 1.0 + EffDofCalc::Get_Calculator().Calc_dlngs_dlnT(T) / 3.0; }

REAL rho_R_at_T(REAL T) { return M_PI * M_PI / 30.0 * f_ge(T) * pow(T, 4); }

REAL T_from_rho_R(REAL rhoR) {
    REAL ge_max = T_gs_ge_Points[N_T_gs_ge_Points - 1][2] +
                  10;                           // using + 10 to further increase ge_max and make sure Tmin < Ttrue
    REAL ge_min = T_gs_ge_Points[0][0] / 10.0;  // using /10 to further decrease ge_min and make sure Tmax > Ttrue
    REAL Tmin = pow(rhoR * 30.0 / M_PI / M_PI / ge_max, 0.25);
    REAL Tmax = pow(rhoR * 30.0 / M_PI / M_PI / ge_min, 0.25);
    REAL tor = 1e-3 * fabs(Tmax - Tmin);
    while (fabs(Tmax - Tmin) > tor) {
        REAL T_test = (Tmax + Tmin) / 2.0;
        REAL rhoR_test = rho_R_at_T(T_test);
        if (rhoR_test < rhoR) {
            Tmin = T_test;
        } else {
            Tmax = T_test;
        }
    }
    return (Tmax - Tmin) / 2.0;
}

}  // namespace EvoEMD
