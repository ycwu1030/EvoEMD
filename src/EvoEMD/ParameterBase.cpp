#include "EvoEMD/ParameterBase.h"

#include <iostream>

namespace EvoEMD {
Parameter_Factory::Parameter_Factory() {}

Parameter_Factory::~Parameter_Factory() {
    Parameter_List::iterator iter_PL = PL.begin();
    for (; iter_PL != PL.end(); iter_PL++) {
        delete iter_PL->second;
    }
}

Parameter_Factory &Parameter_Factory::Get_Parameter_Factory() {
    static Parameter_Factory PF;
    return PF;
}

void Parameter_Factory::Register_Parameter(Parameter_Base *par) {
    std::string name = par->Get_Name();
    Parameter_List::iterator iter = PL.find(name);
    if (iter == PL.end()) {
        PL[name] = par;
    } else {
        std::cout << "Duplicated name: [" << name << "] for parameter" << std::endl;
    }
}

bool Parameter_Factory::Set_Parameter_Value(std::string name, REAL value) {
    Parameter_List::iterator iter = PL.find(name);
    if (iter == PL.end()) {
        return false;
    } else {
        iter->second->Set_Value(value);
        return true;
    }
}

REAL Parameter_Factory::Get_Parameter_Value(std::string name, REAL default_value) {
    Parameter_List::iterator iter_PL = PL.find(name);
    if (iter_PL != PL.end()) {
        return PL[name]->Get_Value();
    }
    std::cout << "Parameter with name: " << name << " is not found" << std::endl;
    return 0;
}

Parameter_Base *Parameter_Factory::Get_Parameter(std::string name) {
    Parameter_List::iterator iter_PL = PL.find(name);
    if (iter_PL == PL.end()) {
        return nullptr;
    }
    return iter_PL->second;
}

}  // namespace EvoEMD
