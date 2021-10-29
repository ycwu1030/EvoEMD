#include "EvoEMD/RealTypes.h"

#include <cmath>

using namespace std;

VD fabs(const VD &input) {
    VD res;
    for (int i = 0; i < input.size(); ++i) {
        res.push_back(std::fabs(input[i]));
    }
    return res;
}

VD operator+(const VD &lhs, const VD &rhs) {
    VD res;
    for (int i = 0; i < lhs.size(); i++) {
        res.push_back(lhs[i] + rhs[i]);
    }
    return res;
}

VD operator+(const VD &lhs, const REAL &cons) {
    VD res;
    for (int i = 0; i < lhs.size(); i++) {
        res.push_back(lhs[i] + cons);
    }
    return res;
}

VD operator+(const REAL &cons, const VD &rhs) {
    VD res;
    for (int i = 0; i < rhs.size(); i++) {
        res.push_back(cons + rhs[i]);
    }
    return res;
}

VD operator-(const VD &lhs, const VD &rhs) {
    VD res;
    for (int i = 0; i < lhs.size(); i++) {
        res.push_back(lhs[i] - rhs[i]);
    }
    return res;
}
VD operator-(const VD &lhs, const REAL &rhs) {
    VD res;
    for (size_t i = 0; i < lhs.size(); i++) {
        res.push_back(lhs[i] - rhs);
    }
    return res;
}
VD operator-(const REAL &lhs, const VD &rhs) {
    VD res;
    for (size_t i = 0; i < rhs.size(); i++) {
        res.push_back(lhs - rhs[i]);
    }
    return res;
}
VD operator-(const VD &rhs) { return 0 - rhs; }

VD operator*(const VD &lhs, const REAL &s) {
    VD res;
    for (int i = 0; i < lhs.size(); i++) {
        res.push_back(lhs[i] * s);
    }
    return res;
}

VD operator*(const REAL &s, const VD &rhs) {
    VD res;
    for (int i = 0; i < rhs.size(); i++) {
        res.push_back(rhs[i] * s);
    }
    return res;
}

REAL operator*(const VD &lhs, const VD &rhs) {
    REAL res = 0;
    for (int i = 0; i < lhs.size(); i++) {
        res += lhs[i] * rhs[i];
    }
    return res;
}
VD operator/(const VD &lhs, const REAL &s) {
    VD res;
    for (int i = 0; i < lhs.size(); ++i) {
        res.push_back(lhs[i] / s);
    }
    return res;
}

VD operator/(const VD &lhs, const VD &rhs) {
    VD res;
    for (int i = 0; i < lhs.size(); ++i) {
        res.push_back(lhs[i] / rhs[i]);
    }
    return res;
}

std::ostream &operator<<(std::ostream &os, const VD &rhs) {
    for (int i = 0; i < rhs.size() - 1; i++) {
        os << rhs[i] << ",";
    }
    os << rhs.back();
    return os;
}
