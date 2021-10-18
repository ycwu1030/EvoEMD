#ifndef __REALTYPES_H__
#define __REALTYPES_H__

#include <iostream>
#include <vector>

#define REAL double
#define VD std::vector<REAL>
#define VVD std::vector<VD>

class ROTATION_NUMBER {
private:
    unsigned value;

public:
    static const unsigned MAX_VERSION_ID = 2048;
    ROTATION_NUMBER() : value(0){};
    ~ROTATION_NUMBER(){};

    ROTATION_NUMBER &operator++();    // ++RN
    ROTATION_NUMBER operator++(int);  // RN++

    bool operator==(const ROTATION_NUMBER &rh);
};

// *
VD fabs(const VD &input);

// * Plus Operator
VD operator+(const VD &lhs, const VD &rhs);
VD operator+(const VD &lhs, const REAL &cons);
VD operator+(const REAL &cons, const VD &rhs);

// * Minus Operator
VD operator-(const VD &lhs, const VD &rhs);
VD operator-(const VD &lhs, const REAL &rhs);
VD operator-(const REAL &lhs, const VD &rhs);
VD operator-(const VD &rhs);

// * Multiply Operator
VD operator*(const VD &lhs, const REAL &s);
VD operator*(const REAL &s, const VD &rhs);
REAL operator*(const VD &lhs, const VD &rhs);  // Scalar Product

// * Divide Operator
VD operator/(const VD &lhs, const VD &rhs);  // elementary-wise divide
VD operator/(const VD &lhs, const REAL &s);

std::ostream &operator<<(std::ostream &os, const VD &rhs);

#endif  //__REALTYPES_H__
