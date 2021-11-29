#include "EvoEMD/HubbleEvolution.h"

#include <algorithm>
#include <cmath>

#include "EvoEMD/Constants.h"
#include "EvoEMD/EffDOF.h"
#include "EvoEMD/HubbleBase.h"
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

Hubble_BE::Hubble_BE() : Parameter_Base("Hubble_BE"), Hubble_Base() {
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
    VD Us = HE.Get_Solution_U();
    VD Y1 = HE.Get_Solution_Y1();
    VD Y2 = HE.Get_Solution_Y2();
    Clean();
    for (int i = 0; i < Us.size(); i++) {
        REAL rhoM = Y1[i] * rhoRr * kappa / pow(Us[i], 3);
        REAL rhoR = Y2[i] * rhoRr * kappa / pow(Us[i], 4);
        REAL Ttemp = T_from_rho_R(rhoR);
        if (i > 0 && Ttemp > List_T.back()) continue;  // Temperature should be decreasing with U increasing.
        REAL HH = sqrt((rhoM + rhoR) / 3.0 / PHY_MP / PHY_MP);
        List_T.push_back(Ttemp);
        List_U.push_back(Us[i]);
        List_H.push_back(HH);
        List_lnT.push_back(log(Ttemp));
        List_lnU.push_back(log(Us[i]));
        List_lnH.push_back(log(HH));
        List_rhoR.push_back(rhoR);
        List_rhoM.push_back(rhoM);
    }
    // Putting all data in T-increasing order
    std::reverse(List_T.begin(), List_T.end());
    std::reverse(List_U.begin(), List_U.end());
    std::reverse(List_H.begin(), List_H.end());
    std::reverse(List_lnT.begin(), List_lnT.end());
    std::reverse(List_lnU.begin(), List_lnU.end());
    std::reverse(List_lnH.begin(), List_lnH.end());
    std::reverse(List_rhoR.begin(), List_rhoR.end());
    std::reverse(List_rhoM.begin(), List_rhoM.end());

    // Interpolation in log scale
    acc_Us = gsl_interp_accel_alloc();
    acc_Hs = gsl_interp_accel_alloc();
    spline_Us = gsl_spline_alloc(gsl_interp_steffen, List_lnT.size());
    spline_Hs = gsl_spline_alloc(gsl_interp_steffen, List_lnT.size());
    gsl_spline_init(spline_Us, List_lnT.data(), List_lnU.data(), List_lnT.size());
    gsl_spline_init(spline_Hs, List_lnT.data(), List_lnH.data(), List_lnT.size());
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
    List_T.clear();
    List_U.clear();
    List_H.clear();
    List_lnT.clear();
    List_lnU.clear();
    List_lnH.clear();
    List_rhoR.clear();
    List_rhoM.clear();
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
    if (T >= List_T.back()) return Get_Hubble_For_RD_at_T(T);
    if (T <= List_T.front()) return Get_Hubble_For_RD_at_T(T);
    return exp(gsl_spline_eval(spline_Hs, log(T), acc_Hs));
}

REAL Hubble_BE::Get_dlna_dlnT_at_T(const REAL T) {
    Get_Value();
    if (T >= List_T.back()) return -1;
    if (T <= List_T.front()) return -1;
    REAL dlnU_dlnT = gsl_spline_eval_deriv(spline_Us, log(T), acc_Us);
    return dlnU_dlnT;
}

void Hubble_BE::Print() { std::cout << "Dummy information for Hubble solved from BE" << std::endl; }

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
