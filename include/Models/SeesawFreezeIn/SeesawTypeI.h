#ifndef _SEESAW_TYPE_I_H_
#define _SEESAW_TYPE_I_H_

#include <Eigen/Dense>
#include <complex>

#include "EvoEMD/ParameterBase.h"

class SeesawTypeI : public EvoEMD::Parameter_Base {
private:
    void Set_PMNS_Matrix();
    void Set_Mass_Ordering();
    void Set_Light_Neutrino_Mass();
    void Set_Heavy_Neutrino_Mass();
    void Set_RHN_Angle();

    void Set_Mixing_Matrix();
    Eigen::Matrix3cd UPMNS_Normal_Order;
    Eigen::Matrix3cd UPMNS_Inverted_Order;
    Eigen::Matrix3cd *UPMNS_Used;

    Eigen::Matrix3cd Mnu_sqrt;
    Eigen::Matrix3cd MNR_sqrt;
    Eigen::Matrix3cd MNR_sqrt_inverse;

    Eigen::Matrix3cd RHN;

    Eigen::Matrix3cd Ynu;
    Eigen::Matrix3cd YdagY;

public:
    SeesawTypeI();
    ~SeesawTypeI(){};

    std::complex<double> Get_Yij(int i, int j);
    std::complex<double> Get_YdagYij(int i, int j);
    virtual void Update_Value(REAL input) override;
};

#endif  //_SEESAW_TYPE_I_H_
