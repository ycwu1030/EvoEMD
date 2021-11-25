#include "EvoEMD/Hubble.h"

#include <cmath>

#include "EvoEMD/Constants.h"
#include "EvoEMD/EffDOF.h"
#include "EvoEMD/ParameterBase.h"

namespace EvoEMD {
Hubble_For_Single_Period::Hubble_For_Single_Period(const REAL Ti, const REAL Tf, const bool Isentropic_in,
                                                   const double beta_R_in)
    : T_start(Ti), T_end(Tf), Isentropic(Isentropic_in), beta_R(beta_R_in) {}

REAL Hubble_For_Single_Period::Get_Hubble_For_RD(const REAL T) {
    REAL geT = f_ge(T);
    return M_PI / 3.0 * sqrt(geT / 10.0) * T * T / PHY_MP;
}

void Hubble_For_Single_Period::Print() const {
    std::cout << "Temperature: [" << T_start << "," << T_end << "], beta_R = " << beta_R
              << ", entropy conservation: " << Isentropic << std::endl;
}

Hubble_RD::Hubble_RD(const REAL T_start, const REAL T_end) : Hubble_For_Single_Period(T_start, T_end) {}

REAL Hubble_RD::Get_Hubble_at_T(const REAL T) { return Get_Hubble_For_RD(T); }

Hubble_EMD::Hubble_EMD(const REAL T_start, const REAL T_end) : Hubble_For_Single_Period(T_start, T_end) {
    HRD_at_T_start = Get_Hubble_For_RD(T_start);
    gs_at_T_start = f_gs(T_start);
}

REAL Hubble_EMD::Get_Hubble_at_T(const REAL T) {
    REAL gs_at_T = f_gs(T);
    return HRD_at_T_start * sqrt(gs_at_T / gs_at_T_start) * pow(T / T_start, 1.5);
}

Hubble_EP::Hubble_EP(const REAL T_start, const REAL T_end)
    : Hubble_For_Single_Period(T_start, T_end, false, 3.0 / 8.0) {
    HRD_at_T_end = Get_Hubble_For_RD(T_end);
    ge_at_T_end = f_ge(T_end);
}

REAL Hubble_EP::Get_Hubble_at_T(const REAL T) {
    REAL ge_at_T = f_ge(T);
    return HRD_at_T_end * ge_at_T / ge_at_T_end * pow(T / T_end, 4);
}

void Hubble_History::Update_Value(REAL input) {
    Ti = RETRIEVE_PARAMETER(Ti)->Get_Value();
    Tr = RETRIEVE_PARAMETER(Tr)->Get_Value();
    TRH = std::max(1e15, Ti * 10);
    Tf = std::min(1e-3, Tr / 10.0);
    Solve_Te();
    for (int i = 0; i < Periods.size(); i++) {
        delete Periods[i];
    }
    Periods.clear();
    Periods.push_back(new Hubble_RD(TRH, Ti));
    Periods.push_back(new Hubble_EMD(Ti, Te));
    Periods.push_back(new Hubble_EP(Te, Tr));
    Periods.push_back(new Hubble_RD(Tr, Tf));
    Temperatures.clear();
    Temperatures.push_back(TRH);
    Temperatures.push_back(Ti);
    Temperatures.push_back(Te);
    Temperatures.push_back(Tr);
    Temperatures.push_back(Tf);
}

Hubble_History::Hubble_History() : Parameter_Base("Hubble") {
    Parameter_Base *ti = RETRIEVE_PARAMETER(Ti);
    Parameter_Base *tr = RETRIEVE_PARAMETER(Tr);
    Register_Dependencies(ti, tr);
    Update_Value(0);
}

Hubble_History::Hubble_History(const Hubble_History &HH)
    : TRH(HH.TRH), Ti(HH.Ti), Te(HH.Te), Tr(HH.Tr), Tf(HH.Tf), Temperatures(HH.Temperatures), Parameter_Base("Hubble") {
    Periods.push_back(new Hubble_RD(TRH, Ti));
    Periods.push_back(new Hubble_EMD(Ti, Te));
    Periods.push_back(new Hubble_EP(Te, Tr));
    Periods.push_back(new Hubble_RD(Tr, Tf));
}

Hubble_History::~Hubble_History() {
    for (int i = 0; i < Periods.size(); i++) {
        delete Periods[i];
    }
}

Hubble_History &Hubble_History::operator=(const Hubble_History &HH) {
    TRH = HH.TRH;
    Ti = HH.Ti;
    Te = HH.Te;
    Tr = HH.Tr;
    Tf = HH.Tf;
    Temperatures = HH.Temperatures;
    for (int i = 0; i < Periods.size(); i++) {
        delete Periods[i];
    }
    Periods.clear();

    Periods.push_back(new Hubble_RD(TRH, Ti));
    Periods.push_back(new Hubble_EMD(Ti, Te));
    Periods.push_back(new Hubble_EP(Te, Tr));
    Periods.push_back(new Hubble_RD(Tr, Tf));

    return *this;
}

void Hubble_History::Solve_Te() {
    // * Te satisfies following function
    // * 5*log(Te) + 2*log(ge(Te)) - log(gs(Te)) == log(Ti) + 4*log(Tr) + log(ge(Tr)) + log(ge(Ti)) - log(gs(Ti));
    // * We solve this equation using binary search, the starting bracket is [Tr,Ti]
    REAL x_max = log(Ti);
    REAL x_min = log(Tr);
    REAL x_eps = 1e-5 * std::fabs(x_max - x_min);
    REAL target_value = log(Ti) + 4 * log(Tr) + log(f_ge(Tr)) + log(f_ge(Ti)) - log(f_gs(Ti));
    REAL test_value;
    REAL x_test;
    while (std::fabs(x_max - x_min) > x_eps) {
        x_test = (x_max + x_min) / 2.0;
        test_value = 5 * x_test + 2 * log(f_ge(exp(x_test))) - log(f_gs(exp(x_test)));
        if (test_value > target_value) {
            x_max = x_test;
        } else {
            x_min = x_test;
        }
    }
    Te = exp((x_max + x_min) / 2.0);
}

int Hubble_History::Get_Period_ID_at_T(const REAL T) {
    if (!updated) {
        Update_Value(0);
    }
    if (T > TRH || T <= Tf) return -1;
    int id_low = 0;
    int id_high = Temperatures.size() - 1;
    int id_mid;
    while (id_high - id_low > 1) {
        id_mid = (id_high + id_low) >> 1;
        if (T <= Temperatures[id_mid]) {
            id_low = id_mid;
        } else {
            id_high = id_mid;
        }
    }
    return id_low;
}

REAL Hubble_History::Get_Hubble_at_T(const REAL T) {
    if (!updated) {
        Update_Value(0);
    }
    return Periods[Get_Period_ID_at_T(T)]->Get_Hubble_at_T(T);
}

double Hubble_History::Get_beta_R_at_T(const REAL T) {
    if (!updated) {
        Update_Value(0);
    }
    return Periods[Get_Period_ID_at_T(T)]->Get_beta_R();
}

Hubble_For_Single_Period *Hubble_History::at(const int pid) {
    if (!updated) {
        Update_Value(0);
    }
    return Periods[pid];
}

Hubble_For_Single_Period *Hubble_History::operator[](const int pid) { return this->at(pid); }

Hubble_History &Hubble_History::Get_Hubble_History() {
    static Hubble_History hist;
    return hist;
}

void Hubble_History::Print_History() {
    if (!updated) {
        Update_Value(0);
    }
    for (int i = 0; i < Periods.size(); i++) {
        Periods[i]->Print();
    }
}

}  // namespace EvoEMD
