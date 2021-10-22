#include "PhysicsParameters.h"

#include <iostream>

Parameter_Factory::Parameter_Factory() {}

Parameter_Factory::~Parameter_Factory() {
    Independent_Parameter_List::iterator iter_IPL = IPL.begin();
    for (; iter_IPL != IPL.end(); iter_IPL++) {
        delete iter_IPL->second;
    }
    Dependent_Parameter_List::iterator iter_DPL = DPL.begin();
    for (; iter_DPL != DPL.end(); iter_DPL++) {
        delete iter_DPL->second;
    }
}

Parameter_Factory &Parameter_Factory::Get_Parameter_Factory() {
    static Parameter_Factory PF;
    return PF;
}

void Parameter_Factory::Register_Parameter(Independent_Parameter *par) {
    std::string name = par->Get_Name();
    Independent_Parameter_List::iterator iter = IPL.find(name);
    if (iter == IPL.end()) {
        IPL[name] = par;
    } else {
        std::cout << "You are duplicated name: [" << name << "] for independent parameter" << std::endl;
    }
}
void Parameter_Factory::Register_Parameter(Dependent_Parameter *par) {
    std::string name = par->Get_Name();
    Dependent_Parameter_List::iterator iter = DPL.find(name);
    if (iter == DPL.end()) {
        DPL[name] = par;
    } else {
        std::cout << "You are duplicated name: [" << name << "] for dependent parameter" << std::endl;
    }
}

bool Parameter_Factory::Set_Independent_Parameter(std::string name, REAL value) {
    Independent_Parameter_List::iterator iter = IPL.find(name);
    if (iter == IPL.end()) {
        return false;
    } else {
        iter->second->Set_Value(value);
        return true;
    }
}

REAL Parameter_Factory::Get_Parameter_Value(std::string name, REAL default_value) {
    Independent_Parameter_List::iterator iter_IPL = IPL.find(name);
    Dependent_Parameter_List::iterator iter_DPL = DPL.find(name);
    if (iter_IPL != IPL.end() && iter_DPL != DPL.end()) {
        std::cout << "You have duplicated parameter name for independent and dependent parameter: " << name
                  << std::endl;
        return IPL[name]->Get_Value();
    }
    if (iter_IPL != IPL.end()) {
        return IPL[name]->Get_Value();
    }
    if (iter_DPL != DPL.end()) {
        return DPL[name]->Get_Value();
    }
    std::cout << "Parameter with name: " << name << " is not found" << std::endl;
    return 0;
}

bool Dependent_Parameter::Is_Updated() {
    bool good;
    Independent_Parameter_Set::iterator iter = IPS.begin();
    for (; iter != IPS.end(); iter++) {
        good = (((iter->second).first)->Get_Version_ID() == (iter->second).second);
        if (!good) return false;
    }
    return true;
}

void Dependent_Parameter::Update_Version_ID() {
    Independent_Parameter_Set::iterator iter = IPS.begin();
    for (; iter != IPS.end(); iter++) {
        (iter->second).second = (((iter->second).first)->Get_Version_ID());
    }
}

REAL Dependent_Parameter::Get_Value() {
    if (!Is_Updated()) {
        Update_Value();
        Update_Version_ID();
    }
    return value;
}
