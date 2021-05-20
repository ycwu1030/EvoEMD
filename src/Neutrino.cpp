#include "Neutrino.h"
#include "spdlog_wrapper.h"
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
using namespace Eigen;

Nu_TypeI_SeeSaw::Nu_TypeI_SeeSaw() : ORDER(NO), UPDATED(false)
{
    SPDLOG_INFO_FILE("Default TypeI SeeSaw Constructor. Normal Ordering");
    Set_PMNS_Matrix();
    Set_Light_Neutrino_Mass();
    Set_Heavy_Neutrino_Mass();
    Set_RHN_Angle();
    Set_Mixing_Matrix();
}

Nu_TypeI_SeeSaw::Nu_TypeI_SeeSaw(MassOrdering od) : ORDER(od), UPDATED(false)
{
    SPDLOG_INFO_FILE("TypeI SeeSaw Constructor. Mass Ordering = {}.", od);
    Set_PMNS_Matrix();
    Set_Light_Neutrino_Mass();
    Set_Heavy_Neutrino_Mass();
    Set_RHN_Angle();
    Set_Mixing_Matrix();
}

void Nu_TypeI_SeeSaw::Set_PMNS_Matrix()
{
    if (ORDER == NO || ORDER == NORMAL_ORDER)
    {
        theta12 = 0.59027;
        theta13 = 0.15027;
        theta23 = 0.86743;
        deltaCP = 3.78736;
        mj_alpha1 = 0.0; // Not available;
        mj_alpha2 = 0.0; // Not available;

        dm221 = 7.39e-5 * eV * eV;
        dm23l = 2.525e-3 * eV * eV;
    }
    else if (ORDER == IO || ORDER == INVERTED_ORDER)
    {
        theta12 = 0.59027;
        theta13 = 0.15097;
        theta23 = 0.86743;
        deltaCP = 4.88692;
        mj_alpha1 = 0.0;
        mj_alpha2 = 0.0;

        dm221 = 7.39e-5 * eV * eV;
        dm23l = -2.512e-3 * eV * eV;
    }
    else
    {
        // ! Un-recognized ordering, using normal order instead
        theta12 = 0.59027;
        theta13 = 0.15027;
        theta23 = 0.86743;
        deltaCP = 3.78736;
        mj_alpha1 = 0.0; // Not available;
        mj_alpha2 = 0.0; // Not available;

        dm221 = 7.39e-5 * eV * eV;
        dm23l = 2.525e-3 * eV * eV;
    }

    complex<double> II = complex<double>(0.0, 1.0);
    complex<double> emicp = cos(deltaCP) - II * sin(deltaCP);
    complex<double> epicp = cos(deltaCP) + II * sin(deltaCP);
    UPMNS(0, 0) = cos(theta12) * cos(theta13);
    UPMNS(0, 1) = cos(theta13) * sin(theta12);
    UPMNS(0, 2) = emicp * sin(theta13);
    UPMNS(1, 0) = -cos(theta23) * sin(theta12) - epicp * cos(theta12) * sin(theta13) * sin(theta23);
    UPMNS(1, 1) = cos(theta12) * cos(theta23) - epicp * sin(theta12) * sin(theta13) * sin(theta23);
    UPMNS(1, 2) = cos(theta13) * sin(theta23);
    UPMNS(2, 0) = sin(theta12) * sin(theta23) - epicp * cos(theta12) * cos(theta23) * sin(theta13);
    UPMNS(2, 1) = -cos(theta12) * sin(theta23) - epicp * cos(theta23) * sin(theta12) * sin(theta13);
    UPMNS(2, 2) = cos(theta13) * cos(theta23);

    UPDATED = false;

    SPDLOG_INFO_FILE("Set Up PMNS Matrix for mass ordering = {}", ORDER);
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(UPMNS(0, 0)), imag(UPMNS(0, 0)), real(UPMNS(0, 1)), imag(UPMNS(0, 1)), real(UPMNS(0, 2)), imag(UPMNS(0, 2)));
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(UPMNS(1, 0)), imag(UPMNS(1, 0)), real(UPMNS(1, 1)), imag(UPMNS(1, 1)), real(UPMNS(1, 2)), imag(UPMNS(1, 2)));
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(UPMNS(2, 0)), imag(UPMNS(2, 0)), real(UPMNS(2, 1)), imag(UPMNS(2, 1)), real(UPMNS(2, 2)), imag(UPMNS(2, 2)));
}

void Nu_TypeI_SeeSaw::Set_Mass_Ordering(MassOrdering od)
{
    ORDER = od;
    Set_PMNS_Matrix();
    UPDATED = false;
}

void Nu_TypeI_SeeSaw::Set_Light_Neutrino_Mass(double m1)
{
    // Light Neutrino:
    if (ORDER == NORMAL_ORDER || ORDER == NO)
    {
        // m1 < m2 < m3
        mnu1 = m1;
        mnu2 = sqrt(mnu1 * mnu1 + dm221);
        mnu3 = sqrt(mnu1 * mnu1 + dm23l);
    }
    else if (ORDER == INVERTED_ORDER || ORDER == IO)
    {
        // m3 < m1 < m2
        mnu3 = m1;
        mnu2 = sqrt(mnu3 * mnu3 - dm23l);
        mnu1 = sqrt(mnu2 * mnu2 - dm221);
    }
    else
    {
        // ! Un-recognized ordering, using normal order instead
        // m1 < m2 < m3
        mnu1 = m1;
        mnu2 = sqrt(mnu1 * mnu1 + dm221);
        mnu3 = sqrt(mnu1 * mnu1 + dm23l);
    }
    Mnu_sqrt = Matrix3cd::Zero();
    Mnu_sqrt(0, 0) = sqrt(mnu1);
    Mnu_sqrt(1, 1) = sqrt(mnu2);
    Mnu_sqrt(2, 2) = sqrt(mnu3);

    UPDATED = false;

    SPDLOG_INFO_FILE("Set up light neutrino mass for mass ordering = {}: m1 = {:+9.8e}, m2 = {:+9.8e}, m3 = {:+9.8e}.", ORDER, mnu1, mnu2, mnu3);
}

void Nu_TypeI_SeeSaw::Set_Heavy_Neutrino_Mass(double m1, double m2, double m3)
{
    vector<double> mm({m1, m2, m3});
    sort(mm.begin(), mm.end());

    MNR1 = mm[0];
    MNR2 = mm[1];
    MNR3 = mm[2];
    MNR_sqrt = Matrix3cd::Zero();
    MNR_sqrt_inverse = Matrix3cd::Zero();
    MNR_sqrt(0, 0) = sqrt(MNR1);
    MNR_sqrt(1, 1) = sqrt(MNR2);
    MNR_sqrt(2, 2) = sqrt(MNR3);
    MNR_sqrt_inverse(0, 0) = 1.0 / sqrt(MNR1);
    MNR_sqrt_inverse(1, 1) = 1.0 / sqrt(MNR2);
    MNR_sqrt_inverse(2, 2) = 1.0 / sqrt(MNR3);
    UPDATED = false;
    SPDLOG_INFO_FILE("Set up heavy neutrino masses: m1 = {:+9.8e}, m2 = {:+9.8e}, m3 = {:+9.8e}.", MNR1, MNR2, MNR3);
}

void Nu_TypeI_SeeSaw::Set_RHN_Angle(double rw12, double iw12, double rw13, double iw13, double rw23, double iw23)
{
    w12_R = rw12;
    w12_I = iw12;
    w13_R = rw13;
    w13_I = iw13;
    w23_R = rw23;
    w23_I = iw23;
    complex<double> II = complex<double>(0.0, 1.0);
    complex<double> c12 = cos(rw12 + II * iw12);
    complex<double> s12 = sin(rw12 + II * iw12);
    complex<double> c13 = cos(rw13 + II * iw13);
    complex<double> s13 = sin(rw13 + II * iw13);
    complex<double> c23 = cos(rw23 + II * iw23);
    complex<double> s23 = sin(rw23 + II * iw23);

    SPDLOG_INFO_FILE("Setting Angles in RHN:");
    SPDLOG_INFO_FILE("w12 = ({:+9.8e}, {:+9.8e})", w12_R, w12_I);
    SPDLOG_INFO_FILE("w13 = ({:+9.8e}, {:+9.8e})", w13_R, w13_I);
    SPDLOG_INFO_FILE("w23 = ({:+9.8e}, {:+9.8e})", w23_R, w23_I);

    RHN(0, 0) = c12 * c13;
    RHN(0, 1) = s12 * c13;
    RHN(0, 2) = s13;
    RHN(1, 0) = -c23 * s12 - c12 * s13 * s23;
    RHN(1, 1) = c12 * c23 - s12 * s13 * s23;
    RHN(1, 2) = c13 * s23;
    RHN(2, 0) = s12 * s23 - c12 * c23 * s13;
    RHN(2, 1) = -c12 * s23 - c23 * s12 * s13;
    RHN(2, 2) = c13 * c23;

    SPDLOG_INFO_FILE("Set Up R Matrix:");
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(RHN(0, 0)), imag(RHN(0, 0)), real(RHN(0, 1)), imag(RHN(0, 1)), real(RHN(0, 2)), imag(RHN(0, 2)));
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(RHN(1, 0)), imag(RHN(1, 0)), real(RHN(1, 1)), imag(RHN(1, 1)), real(RHN(1, 2)), imag(RHN(1, 2)));
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(RHN(2, 0)), imag(RHN(2, 0)), real(RHN(2, 1)), imag(RHN(2, 1)), real(RHN(2, 2)), imag(RHN(2, 2)));

    UPDATED = false;
}

void Nu_TypeI_SeeSaw::Set_Mixing_Matrix()
{
    complex<double> II = complex<double>(0.0, 1.0);
    Ynu = II * sqrt(2.0) / 246.0 * UPMNS * Mnu_sqrt * RHN * MNR_sqrt;
    YdagY = Ynu.adjoint() * Ynu;
    SPDLOG_INFO_FILE("Updating the Yukawa Coupling matrix:");
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(Ynu(0, 0)), imag(Ynu(0, 0)), real(Ynu(0, 1)), imag(Ynu(0, 1)), real(Ynu(0, 2)), imag(Ynu(0, 2)));
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(Ynu(1, 0)), imag(Ynu(1, 0)), real(Ynu(1, 1)), imag(Ynu(1, 1)), real(Ynu(1, 2)), imag(Ynu(1, 2)));
    SPDLOG_INFO_FILE("({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})  ({:+9.8e},{:+9.8e})", real(Ynu(2, 0)), imag(Ynu(2, 0)), real(Ynu(2, 1)), imag(Ynu(2, 1)), real(Ynu(2, 2)), imag(Ynu(2, 2)));
    UPDATED = true;
}

complex<double> Nu_TypeI_SeeSaw::Get_Yij(int i, int j)
{
    if (!UPDATED)
        Set_Mixing_Matrix();
    return Ynu(i, j);
}
complex<double> Nu_TypeI_SeeSaw::Get_Yij(int i, int j)
{
    if (!UPDATED)
        Set_Mixing_Matrix();
    return YdagY(i, j);
}

complex<double> Nu_TypeI_SeeSaw::Get_UPMNSij(int i, int j)
{
    return UPMNS(i, j);
}
