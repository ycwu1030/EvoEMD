#ifndef __NEUTRINO_H__
#define __NEUTRINO_H__

#include <Eigen/Dense>
#include <complex>

#include "Particles.h"
#include "Physics_Constants.h"

class Nu_TypeI_SeeSaw : public Particle_Client {
public:
    typedef enum {
        NORMAL_ORDER = 0,
        INVERTED_ORDER = 1,
        NO = 2,  // * Alternative name for normal order
        IO = 3   // * Alternative name for inverted order
    } MassOrdering;
    Nu_TypeI_SeeSaw();
    // Nu_TypeI_SeeSaw(MassOrdering od);
    ~Nu_TypeI_SeeSaw() = default;
    static Nu_TypeI_SeeSaw& Get_Neutrino_Model();

    virtual void Update_Particle_Info() override;

    /*
     * @brief set the mass ordering and then setup the PMNS matrix.
     */
    void Set_Mass_Ordering(MassOrdering od);

    /*
     * @brief set the light neutrino mass
     *        Only the lightest one is needed.
     *        The other two are determined by the mass ording and the mass difference measurement
     * @param m1, the lightest neutrino mass
     */
    void Set_Light_Neutrino_Mass(double m1 = 0.1 * eV);

    /*
     * @brief set the heavy neutrino mass
     *        m1 <= m2 <= m3
     * @param m1, mass for one of the heavy neutrino
     * @param m2, mass for another heavy neutrino
     * @param m3, mass for the last heavy neutrino
     */
    void Set_Heavy_Neutrino_Mass(double m1 = 3000 * GeV, double m2 = 3000 * GeV, double m3 = 3000 * GeV);

    /*
     * @brief set the angle in the matrix R
     * @param rw12,iw12, the real/imaginary part of w12
     * @param rw13,iw13, the real/imaginary part of w13
     * @param rw23,iw23, the real/imaginary part of w23
     */
    void Set_RHN_Angle(double rw12 = 0, double iw12 = 0, double rw13 = 0, double iw13 = 0, double rw23 = 0,
                       double iw23 = 0);

    /*
     * @brief get the (i,j) element of the Yukawa Coupling Matrix
     */
    std::complex<double> Get_Yij(int i, int j);
    std::complex<double> Get_YdagYij(int i, int j);

    /*
     * @brief get the (i,j) element of the Yukawa Coupling Matrix
     */
    std::complex<double> Get_UPMNSij(int i, int j);

    double Get_NR_Mass(int i);

private:
    // * Flag
    bool UPDATED;
    // * Angles in PMNS
    double theta12;
    double theta13;
    double theta23;
    double deltaCP;

    // * Possible Majorana phases
    double mj_alpha1;
    double mj_alpha2;

    // * PMNS Matrix
    Eigen::Matrix3cd UPMNS;

    // * Angles in R matrix,
    double w12_R, w12_I;
    double w13_R, w13_I;
    double w23_R, w23_I;

    // * Complex Orthogonal Matrix
    Eigen::Matrix3cd RHN;

    // * Light neutrino mass square difference;
    MassOrdering ORDER;
    double dm221;  // * m2^2 - m1^2
    double dm23l;  // * m3^2 - ml^2; For NO, l=1 and dm23l>0; For IO, l=2 and dm23l<0
    double mnu1;
    double mnu2;
    double mnu3;

    // * Diagonal sqrt(Mass) matrix for light neutrino
    Eigen::Matrix3cd Mnu_sqrt;

    // * Heavy Neutrino masses:
    double MNR1;
    double MNR2;
    double MNR3;

    // * Diagonal sqrt(Mass) matrix for heavy neutrino
    Eigen::Matrix3cd MNR_sqrt;
    Eigen::Matrix3cd MNR_sqrt_inverse;

    // * The Yukawa coupling matrix  Ynu L.H.nuR
    Eigen::Matrix3cd Ynu;    // * = i sqrt2/vev UPMNS.Mnu_sqrt.RHN.MNR_sqrt;
    Eigen::Matrix3cd YdagY;  // * Y^dagger*Y

    /*
     * @brief: Set up the PMNS matrix based on current mass ordering
     *         Will also set the mass square difference
     */
    void Set_PMNS_Matrix();

    /*
     * @brief: Set up the matrix Ynu
     *         This function will be called based on flag: UPDATED
     */
    void Set_Mixing_Matrix();
};

#endif  //__NEUTRINO_H__
