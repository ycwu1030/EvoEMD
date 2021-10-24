#ifndef _PARAMETER_BASE_H_
#define _PARAMETER_BASE_H_

#include <map>
#include <set>
#include <string>

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

class Parameter_Base {
protected:
    std::string name;
    REAL value;
    bool updated;
    std::set<Parameter_Base *> underlying_parameters;  // Parameter in this set depends on current parameter

public:
    Parameter_Base(std::string par_name) : name(par_name), value(0), updated(true){};
    virtual ~Parameter_Base(){};

    bool Is_Independent() { return (underlying_parameters.size() == 0); }
    void Claim_Dependence(Parameter_Base *par) { underlying_parameters.insert(par); }
    void Notify() {
        updated = false;
        for (auto &&bp : underlying_parameters) {
            bp->Notify();
        }
    }

    // * Any parameter derived from this Base_Parameter should re-implement this function
    virtual void Update_Value(REAL input) = 0;

    void Set_Value(REAL input = 0) {
        Notify();
        Update_Value(input);
        updated = true;
    }
    REAL Get_Value() {
        if (updated) {
            return value;
        } else {
            Set_Value();
            return value;
        }
    }
    std::string Get_Name() const { return name; }
};

class Free_Parameter : public Parameter_Base {
public:
    Free_Parameter(std::string name, REAL default_value = 0) : Parameter_Base(name) {
        value = default_value;
        updated = true;
    }
    ~Free_Parameter(){};

    virtual void Update_Value(REAL input) override { value = input; }
};

class Parameter_Factory {
private:
    Parameter_Factory();
    ~Parameter_Factory();
    typedef std::map<std::string, Parameter_Base *> Parameter_List;
    Parameter_List PL;

public:
    static Parameter_Factory &Get_Parameter_Factory();

    void Register_Parameter(Parameter_Base *par);
    bool Set_Parameter_Value(std::string name, REAL value);
    REAL Get_Parameter_Value(std::string name, REAL default_value = 0);
    Parameter_Base *Get_Parameter(std::string name);
    void List_Parameters();
};
}  // namespace EvoEMD

class Register_Parameter {
public:
    Register_Parameter(EvoEMD::Parameter_Base *par) {
        EvoEMD::Parameter_Factory::Get_Parameter_Factory().Register_Parameter(par);
    };
};

#define REGISTER_PARAMETER(paramName) Register_Parameter g_register_parameter_##paramName(new paramName)
#define DECLARE_FREE_PARAMETER(paramName, value)                   \
    class param_##paramName : public Free_Parameter {              \
    public:                                                        \
        param_##paramName() : Free_Parameter(#paramName, value){}; \
    };                                                             \
    REGISTER_PARAMETER(param_##paramName)

#define RETRIVE_PARAMETER(paramName) EvoEMD::Parameter_Factory::Get_Parameter_Factory().Get_Parameter(#paramName)

#endif  //_PARAMETER_BASE_H_
