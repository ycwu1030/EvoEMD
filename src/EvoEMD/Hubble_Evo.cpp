#include "EvoEMD/Hubble_Evo.h"

#include <algorithm>
#include <cmath>

#include "EvoEMD/Constants.h"
#include "EvoEMD/EffDOF.h"
#include "gsl/gsl_spline.h"

namespace EvoEMD {

Hubble_Evolution::Hubble_Evolution() : rk(this), BR(1) { Set_DOF(2); }

void Hubble_Evolution::Solve(REAL Uini, REAL Y1ini, REAL Y2ini, REAL BR_in) {
    BR = BR_in;
    X_BEGIN = log(Uini);
    X_END = X_BEGIN + 50;
    Y_BEGIN[0] = Y1ini;
    Y_BEGIN[1] = Y2ini;
    rk.Solve(1e-3, 1e-4);
    X = rk.Get_Solution_X();
    VVD Y = rk.Get_Solution_Y();
    U = X;
    Y1.clear();
    Y2.clear();
    for (int i = 0; i < U.size(); i++) {
        U[i] = exp(X[i]);
        Y1.push_back(Y[i][0]);
        Y2.push_back(Y[i][1]);
    }
}

Hubble_BE::Hubble_BE() : Parameter_Base("Hubble_BE"), Hubble_For_Single_Period(-1, -1, false, 1) {
    Parameter_Base *ptr_ti = RETRIEVE_PARAMETER(Ti);
    Parameter_Base *ptr_tr = RETRIEVE_PARAMETER(Tr);
    Parameter_Base *ptr_br = RETRIEVE_PARAMETER(BR);

    Register_Dependencies(ptr_ti, ptr_tr, ptr_br);
    acc_Us = nullptr;
    acc_Hs = nullptr;
    spline_Us = nullptr;
    spline_Hs = nullptr;
}

void Hubble_BE::Update_Value(REAL input) {
    REAL ti = RETRIEVE_PARAMETER(Ti)->Get_Value();
    REAL tr = RETRIEVE_PARAMETER(Tr)->Get_Value();
    REAL br = RETRIEVE_PARAMETER(BR)->Get_Value();
    REAL rhoRi = M_PI * M_PI / 30.0 * f_ge(ti) * pow(ti, 4);
    REAL rhoRr = M_PI * M_PI / 30.0 * f_ge(tr) * pow(tr, 4);
    REAL tr_obtained = 0;
    REAL tor = 1e-3;
    REAL U0 = 1.0;
    REAL log_kappa_max = log10(rhoRi / rhoRr);
    REAL log_kappa_min = -4;
    while (fabs(log_kappa_max - log_kappa_min) > tor) {
        REAL log_kappa = (log_kappa_max + log_kappa_min) / 2.0;
        REAL kappa = pow(10, log_kappa);
        REAL Y1init = rhoRi / rhoRr * pow(U0, 3) / kappa;
        REAL Y2init = rhoRi / rhoRr * pow(U0, 4) / kappa;
        HE.Solve(U0, Y1init, Y2init, br);
        tr_obtained = Calc_Tr(HE.Get_Solution_U(), HE.Get_Solution_Y1(), HE.Get_Solution_Y2(), kappa * rhoRr);
        if (tr_obtained < tr) {
            log_kappa_min = log_kappa;
        } else {
            log_kappa_max = log_kappa;
        }
    }
    REAL log_kappa = (log_kappa_max + log_kappa_min) / 2.0;
    REAL kappa = pow(10, log_kappa);
    REAL Y1init = rhoRi / rhoRr * pow(U0, 3) / kappa;
    REAL Y2init = rhoRi / rhoRr * pow(U0, 4) / kappa;
    HE.Solve(U0 / 100.0, Y1init, Y2init, br);
    List_U = HE.Get_Solution_U();
    VD Y1 = HE.Get_Solution_Y1();
    VD Y2 = HE.Get_Solution_Y2();
    List_T.clear();
    List_H.clear();
    List_rhoR.clear();
    List_rhoM.clear();
    for (int i = 0; i < List_U.size(); i++) {
        double rhoM = Y1[i] * rhoRr * kappa / pow(List_U[i], 3);
        double rhoR = Y2[i] * rhoRr * kappa / pow(List_U[i], 4);
        List_T.push_back(T_from_rho_R(rhoR));
        List_H.push_back(sqrt((rhoM + rhoR) / 3.0 / PHY_MP / PHY_MP));
        List_rhoR.push_back(rhoR);
        List_rhoM.push_back(rhoM);
    }
    Clean();
    acc_Us = gsl_interp_accel_alloc();
    acc_Hs = gsl_interp_accel_alloc();
    spline_Us = gsl_spline_alloc(gsl_interp_steffen, List_T.size());
    spline_Hs = gsl_spline_alloc(gsl_interp_steffen, List_T.size());
    gsl_spline_init(spline_Us, List_T.data(), List_U.data(), List_T.size());
    gsl_spline_init(spline_Hs, List_T.data(), List_H.data(), List_T.size());
}

REAL Hubble_BE::Calc_Tr(VD U, VD Y1, VD Y2, REAL krhoRr) {
    VD ratio;
    for (int i = 0; i < U.size(); i++) {
        ratio.push_back(Y1[i] / Y2[i] * U[i]);
    }
    auto iter = std::find_if(ratio.rbegin(), ratio.rend(), [](REAL x) { return x >= 1.0; });
    int dist = std::distance(ratio.rbegin(), iter);
    auto iterY2 = Y2.rbegin() + dist;
    auto iterU = U.rbegin() + dist;
    REAL tau = ((*iter) - 1) / ((*iter) - (*(iter - 1)));
    REAL Y2high = *iterY2;
    REAL Y2low = *(iterY2 - 1);
    REAL Y2res = Y2high - tau * (Y2high - Y2low);
    REAL Uhigh = *iterU;
    REAL Ulow = *(iterU - 1);
    REAL Ures = Uhigh - tau * (Uhigh - Ulow);
    REAL rhoRres = Y2res * krhoRr / pow(Ures, 4);
    return T_from_rho_R(rhoRres);
}

void Hubble_BE::Clean() {
    if (spline_Us) {
        gsl_spline_free(spline_Us);
    }
    if (spline_Hs) {
        gsl_spline_free(spline_Hs);
    }
    if (acc_Us) {
        gsl_interp_accel_free(acc_Us);
    }
    if (acc_Hs) {
        gsl_interp_accel_free(acc_Hs);
    }
}

REAL Hubble_BE::Get_Hubble_at_T(const REAL T) {
    Get_Value();
    if (T >= List_T.front()) return Get_Hubble_For_RD(T);
    if (T <= List_T.back()) return Get_Hubble_For_RD(T);
    return gsl_spline_eval(spline_Hs, T, acc_Hs);
}

REAL Hubble_BE::Get_dlna_dlnT_at_T(const REAL T) {
    Get_Value();
    if (T >= List_T.front()) return -1;
    if (T <= List_T.back()) return -1;
    REAL dU_dT = gsl_spline_eval_deriv(spline_Us, T, acc_Us);
    REAL U = gsl_spline_eval(spline_Us, T, acc_Us);
    return dU_dT / U * T;
}

Hubble_BE &Hubble_BE::Get_Hubble_Calculator() {
    static Hubble_BE HBE;
    return HBE;
}

VD Hubble_BE::Get_T_List() {
    Get_Value();
    return List_T;
}

VD Hubble_BE::Get_U_List() {
    Get_Value();
    return List_U;
}

VD Hubble_BE::Get_H_List() {
    Get_Value();
    return List_H;
}

VD Hubble_BE::Get_rhoR_List() {
    Get_Value();
    return List_rhoR;
}

VD Hubble_BE::Get_rhoM_List() {
    Get_Value();
    return List_rhoM;
}

}  // namespace EvoEMD
