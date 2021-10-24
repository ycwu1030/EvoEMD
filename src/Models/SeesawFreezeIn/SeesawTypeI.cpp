#include "Models/SeesawFreezeIn/SeesawTypeI.h"

#include <cmath>

#include "EvoEMD/Constants.h"

void SeesawTypeI::Set_PMNS_Matrix() {
    // For Normal Order
    double theta12 = 0.59027;
    double theta13 = 0.15027;
    double theta23 = 0.86743;
    double deltaCP = 3.78736;
    double mj_alpha1 = 0.0;  // Not available;
    double mj_alpha2 = 0.0;  // Not available;

    std::complex<double> II = std::complex<double>(0.0, 1.0);
    std::complex<double> emicp = cos(deltaCP) - II * sin(deltaCP);
    std::complex<double> epicp = cos(deltaCP) + II * sin(deltaCP);
    UPMNS_Normal_Order(0, 0) = cos(theta12) * cos(theta13);
    UPMNS_Normal_Order(0, 1) = cos(theta13) * sin(theta12);
    UPMNS_Normal_Order(0, 2) = emicp * sin(theta13);
    UPMNS_Normal_Order(1, 0) = -cos(theta23) * sin(theta12) - epicp * cos(theta12) * sin(theta13) * sin(theta23);
    UPMNS_Normal_Order(1, 1) = cos(theta12) * cos(theta23) - epicp * sin(theta12) * sin(theta13) * sin(theta23);
    UPMNS_Normal_Order(1, 2) = cos(theta13) * sin(theta23);
    UPMNS_Normal_Order(2, 0) = sin(theta12) * sin(theta23) - epicp * cos(theta12) * cos(theta23) * sin(theta13);
    UPMNS_Normal_Order(2, 1) = -cos(theta12) * sin(theta23) - epicp * cos(theta23) * sin(theta12) * sin(theta13);
    UPMNS_Normal_Order(2, 2) = cos(theta13) * cos(theta23);

    // For inverted order
    theta12 = 0.59027;
    theta13 = 0.15097;
    theta23 = 0.86743;
    deltaCP = 4.88692;
    mj_alpha1 = 0.0;
    mj_alpha2 = 0.0;

    emicp = cos(deltaCP) - II * sin(deltaCP);
    epicp = cos(deltaCP) + II * sin(deltaCP);
    UPMNS_Inverted_Order(0, 0) = cos(theta12) * cos(theta13);
    UPMNS_Inverted_Order(0, 1) = cos(theta13) * sin(theta12);
    UPMNS_Inverted_Order(0, 2) = emicp * sin(theta13);
    UPMNS_Inverted_Order(1, 0) = -cos(theta23) * sin(theta12) - epicp * cos(theta12) * sin(theta13) * sin(theta23);
    UPMNS_Inverted_Order(1, 1) = cos(theta12) * cos(theta23) - epicp * sin(theta12) * sin(theta13) * sin(theta23);
    UPMNS_Inverted_Order(1, 2) = cos(theta13) * sin(theta23);
    UPMNS_Inverted_Order(2, 0) = sin(theta12) * sin(theta23) - epicp * cos(theta12) * cos(theta23) * sin(theta13);
    UPMNS_Inverted_Order(2, 1) = -cos(theta12) * sin(theta23) - epicp * cos(theta23) * sin(theta12) * sin(theta13);
    UPMNS_Inverted_Order(2, 2) = cos(theta13) * cos(theta23);
}
void SeesawTypeI::Set_Mass_Ordering() {
    Parameter_Base *od = RETRIVE_PARAMETER(nu_order);
    REAL vod = od->Get_Value();
    if (vod >= 0) {
        UPMNS_Used = &UPMNS_Normal_Order;
    } else {
        UPMNS_Used = &UPMNS_Inverted_Order;
    }
}
void SeesawTypeI::Set_Light_Neutrino_Mass() {
    Parameter_Base *od = RETRIVE_PARAMETER(nu_order);
    Parameter_Base *lm = RETRIVE_PARAMETER(Mnu1);
    REAL vod = od->Get_Value();
    REAL mnu1;
    REAL mnu2;
    REAL mnu3;
    static const REAL dm221_no = 7.39e-5 * eV * eV;
    static const REAL dm23l_no = 2.525e-3 * eV * eV;
    static const REAL dm221_io = 7.39e-5 * eV * eV;
    static const REAL dm23l_io = -2.512e-3 * eV * eV;
    if (vod >= 0) {
        mnu1 = lm->Get_Value();
        mnu2 = sqrt(mnu1 * mnu1 + dm221_no);
        mnu3 = sqrt(mnu1 * mnu1 + dm23l_no);
    } else {
        mnu3 = lm->Get_Value();
        mnu2 = sqrt(mnu3 * mnu3 - dm23l_io);
        mnu1 = sqrt(mnu2 * mnu2 - dm221_io);
    }
    Mnu_sqrt.setZero();
    Mnu_sqrt(0, 0) = sqrt(mnu1);
    Mnu_sqrt(1, 1) = sqrt(mnu2);
    Mnu_sqrt(2, 2) = sqrt(mnu3);
}
void SeesawTypeI::Set_Heavy_Neutrino_Mass() {
    REAL mn1 = RETRIVE_PARAMETER(MN1)->Get_Value();
    REAL mn2 = RETRIVE_PARAMETER(MN2)->Get_Value();
    REAL mn3 = RETRIVE_PARAMETER(MN3)->Get_Value();

    MNR_sqrt.setZero();
    MNR_sqrt_inverse.setZero();
    MNR_sqrt(0, 0) = sqrt(mn1);
    MNR_sqrt(1, 1) = sqrt(mn2);
    MNR_sqrt(2, 2) = sqrt(mn3);
    MNR_sqrt_inverse(0, 0) = 1.0 / sqrt(mn1);
    MNR_sqrt_inverse(1, 1) = 1.0 / sqrt(mn2);
    MNR_sqrt_inverse(2, 2) = 1.0 / sqrt(mn3);
}

void SeesawTypeI::Set_RHN_Angle() {
    REAL w12R = RETRIVE_PARAMETER(rw12)->Get_Value();
    REAL w12I = RETRIVE_PARAMETER(iw12)->Get_Value();
    REAL w13R = RETRIVE_PARAMETER(rw13)->Get_Value();
    REAL w13I = RETRIVE_PARAMETER(iw13)->Get_Value();
    REAL w23R = RETRIVE_PARAMETER(rw23)->Get_Value();
    REAL w23I = RETRIVE_PARAMETER(iw23)->Get_Value();

    std::complex<double> II = std::complex<double>(0.0, 1.0);
    std::complex<double> c12 = cos(w12R + II * w12I);
    std::complex<double> s12 = sin(w12R + II * w12I);
    std::complex<double> c13 = cos(w13R + II * w13I);
    std::complex<double> s13 = sin(w13R + II * w13I);
    std::complex<double> c23 = cos(w23R + II * w23I);
    std::complex<double> s23 = sin(w23R + II * w23I);

    RHN(0, 0) = c12 * c13;
    RHN(0, 1) = s12 * c13;
    RHN(0, 2) = s13;
    RHN(1, 0) = -c23 * s12 - c12 * s13 * s23;
    RHN(1, 1) = c12 * c23 - s12 * s13 * s23;
    RHN(1, 2) = c13 * s23;
    RHN(2, 0) = s12 * s23 - c12 * c23 * s13;
    RHN(2, 1) = -c12 * s23 - c23 * s12 * s13;
    RHN(2, 2) = c13 * c23;
}

void SeesawTypeI::Set_Mixing_Matrix() {
    std::complex<double> II = std::complex<double>(0.0, 1.0);
    Ynu = II * sqrt(2.0) / 246.0 * (*UPMNS_Used) * Mnu_sqrt * RHN * MNR_sqrt;
    YdagY = Ynu.adjoint() * Ynu;
}

SeesawTypeI::SeesawTypeI() : EvoEMD::Parameter_Base("Seesaw_Parameters") {
    Set_PMNS_Matrix();
    RETRIVE_PARAMETER(MN1)->Claim_Dependence(this);
    RETRIVE_PARAMETER(MN2)->Claim_Dependence(this);
    RETRIVE_PARAMETER(MN3)->Claim_Dependence(this);
    RETRIVE_PARAMETER(Mnu1)->Claim_Dependence(this);
    RETRIVE_PARAMETER(rw12)->Claim_Dependence(this);
    RETRIVE_PARAMETER(iw12)->Claim_Dependence(this);
    RETRIVE_PARAMETER(rw23)->Claim_Dependence(this);
    RETRIVE_PARAMETER(iw23)->Claim_Dependence(this);
    RETRIVE_PARAMETER(rw13)->Claim_Dependence(this);
    RETRIVE_PARAMETER(iw13)->Claim_Dependence(this);
    RETRIVE_PARAMETER(nu_order)->Claim_Dependence(this);
    Set_PMNS_Matrix();
    Update_Value(0);
}

void SeesawTypeI::Update_Value(REAL input) {
    value = input;  // Dummy value, avoid un-used argument
    Set_Mass_Ordering();
    Set_Light_Neutrino_Mass();
    Set_Heavy_Neutrino_Mass();
    Set_RHN_Angle();
    Set_Mixing_Matrix();
}

std::complex<double> SeesawTypeI::Get_Yij(int i, int j) {
    if (updated) {
        return Ynu(i, j);
    }
    Set_Value();
    return Ynu(i, j);
}

std::complex<double> SeesawTypeI::Get_YdagYij(int i, int j) {
    if (updated) {
        return YdagY(i, j);
    }
    Set_Value();
    return YdagY(i, j);
}
