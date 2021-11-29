#include <cmath>

#include "EvoEMD/Constants.h"
#include "EvoEMD/EffDOF.h"
#include "EvoEMD/HubbleBase.h"
#include "EvoEMD/ParameterBase.h"

namespace EvoEMD {

REAL Get_Hubble_For_RD_at_T(const REAL T) {
    REAL geT = f_ge(T);
    return M_PI / 3.0 * sqrt(geT / 10.0) * T * T / PHY_MP;
}

Hubble_Splitting_Period::Hubble_Splitting_Period(const REAL T_begin, const REAL T_end, const double beta_R_in)
    : T_BEGIN(T_begin), T_END(T_end), beta_R(beta_R_in), Hubble_Base() {}

REAL Hubble_Splitting_Period::Get_Hubble_at_T(const REAL T) { return Get_Hubble_For_RD_at_T(T); }

REAL Hubble_Splitting_Period::Get_dlna_dlnT_at_T(const REAL T) { return -f_ge_star(T) / beta_R; }

void Hubble_Splitting_Period::Print() {
    std::cout << "Temperature: [" << T_BEGIN << "," << T_END << "], beta_R: " << beta_R << std::endl;
}

Hubble_RD::Hubble_RD(const REAL T_begin, const REAL T_end) : Hubble_Splitting_Period(T_begin, T_end) {}

Hubble_EMD::Hubble_EMD(const REAL T_begin, const REAL T_end) : Hubble_Splitting_Period(T_begin, T_end) {
    HRD_at_T_Begin = Get_Hubble_For_RD_at_T(T_begin);
    gs_at_T_Begin = f_gs(T_begin);
}

REAL Hubble_EMD::Get_Hubble_at_T(const REAL T) {
    REAL gs_at_T = f_gs(T);
    return HRD_at_T_Begin * sqrt(gs_at_T / gs_at_T_Begin) * pow(T / T_BEGIN, 1.5);
}

Hubble_EP::Hubble_EP(const REAL T_begin, const REAL T_end) : Hubble_Splitting_Period(T_begin, T_end, 3.0 / 8.0) {
    HRD_at_T_end = Get_Hubble_For_RD_at_T(T_end);
    ge_at_T_end = f_ge(T_end);
}

REAL Hubble_EP::Get_Hubble_at_T(const REAL T) {
    REAL ge_at_T = f_ge(T);
    return HRD_at_T_end * ge_at_T / ge_at_T_end * pow(T / T_END, 4);
}

void Hubble_Splitting::Update_Value(REAL input) {
    Ti = RETRIEVE_PARAMETER(Ti)->Get_Value();
    Tr = RETRIEVE_PARAMETER(Tr)->Get_Value();
    Solve_Te();
    for (int i = 0; i < Periods.size(); i++) {
        delete Periods[i];
    }
    Periods.clear();
    Periods.push_back(new Hubble_RD(-1, Ti));
    Periods.push_back(new Hubble_EMD(Ti, Te));
    Periods.push_back(new Hubble_EP(Te, Tr));
    Periods.push_back(new Hubble_RD(Tr, -1));
    Temperatures.clear();
    Temperatures.push_back(Ti);
    Temperatures.push_back(Te);
    Temperatures.push_back(Tr);
}

Hubble_Splitting::Hubble_Splitting() : Parameter_Base("Hubble_Splitting") {
    Parameter_Base *ti = RETRIEVE_PARAMETER(Ti);
    Parameter_Base *tr = RETRIEVE_PARAMETER(Tr);
    Register_Dependencies(ti, tr);
    Update_Value(0);
}

Hubble_Splitting::Hubble_Splitting(const Hubble_Splitting &HH)
    : TRH(HH.TRH), Ti(HH.Ti), Te(HH.Te), Tr(HH.Tr), Tf(HH.Tf), Temperatures(HH.Temperatures), Parameter_Base("Hubble") {
    Periods.push_back(new Hubble_RD(-1, Ti));
    Periods.push_back(new Hubble_EMD(Ti, Te));
    Periods.push_back(new Hubble_EP(Te, Tr));
    Periods.push_back(new Hubble_RD(Tr, -1));
}

Hubble_Splitting::~Hubble_Splitting() {
    for (int i = 0; i < Periods.size(); i++) {
        delete Periods[i];
    }
}

Hubble_Splitting &Hubble_Splitting::operator=(const Hubble_Splitting &HH) {
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

    Periods.push_back(new Hubble_RD(-1, Ti));
    Periods.push_back(new Hubble_EMD(Ti, Te));
    Periods.push_back(new Hubble_EP(Te, Tr));
    Periods.push_back(new Hubble_RD(Tr, -1));

    return *this;
}

void Hubble_Splitting::Solve_Te() {
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

int Hubble_Splitting::Get_Period_ID_at_T(const REAL T) {
    if (!updated) {
        Update_Value(0);
    }
    if (T >= Ti) return 0;  // ERD
    if (T <= Tr) return 3;  // RD
    if (T > Te) return 1;   // EMD
    return 2;               // EP
    // int id_low = 0;
    // int id_high = Temperatures.size() - 1;
    // int id_mid;
    // while (id_high - id_low > 1) {
    //     id_mid = (id_high + id_low) >> 1;
    //     if (T <= Temperatures[id_mid]) {
    //         id_low = id_mid;
    //     } else {
    //         id_high = id_mid;
    //     }
    // }
    // return id_low;
}

REAL Hubble_Splitting::Get_Hubble_at_T(const REAL T) {
    if (!updated) {
        Update_Value(0);
    }
    return Periods[Get_Period_ID_at_T(T)]->Get_Hubble_at_T(T);
}

REAL Hubble_Splitting::Get_dlna_dlnT_at_T(const REAL T) {
    if (!updated) {
        Update_Value(0);
    }
    return Periods[Get_Period_ID_at_T(T)]->Get_dlna_dlnT_at_T(T);
}

Hubble_Base *Hubble_Splitting::at(const int pid) {
    if (!updated) {
        Update_Value(0);
    }
    return Periods[pid];
}

Hubble_Base *Hubble_Splitting::operator[](const int pid) { return this->at(pid); }

void Hubble_Splitting::Print() {
    if (!updated) {
        Update_Value(0);
    }
    for (int i = 0; i < Periods.size(); i++) {
        Periods[i]->Print();
    }
}

Hubble_Factory::Hubble_Factory() {
    HCL[0] = new Hubble_Splitting();
    HCL[1] = new Hubble_BE();
}

Hubble_Factory::~Hubble_Factory() {
    for (auto &&hh : HCL) {
        delete hh.second;
    }
}

Hubble_Factory &Hubble_Factory::Get_Hubble_Factory() {
    static Hubble_Factory HF;
    return HF;
}

Hubble_Base *Hubble_Factory::Get_Hubble_Calculator(int id) {
    unsigned int picked_id;
    if (id < 0) {
        id = RETRIEVE_PARAMETER(HubbleID)->Get_Value();
    }
    picked_id = id;
    auto iter = HCL.find(picked_id);
    if (iter == HCL.end()) iter = HCL.begin();
    // iter->second->Print();
    return iter->second;
}

void Hubble_Factory::Register_Hubble_Calculator(int id, Hubble_Base *ptr) {
    if (HCL.count(id) == 1) {
        std::cout << "Duplicated ID for Hubble Calculator encountered. The calculator will be replaced" << std::endl;
    }
    HCL[id] = ptr;
}

}  // namespace EvoEMD
