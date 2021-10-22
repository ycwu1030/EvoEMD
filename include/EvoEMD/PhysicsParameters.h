#ifndef _PHYSICS_PARAMETERS_H_
#define _PHYSICS_PARAMETERS_H_

#include <map>
#include <set>
#include <string>

#include "EvoEMD/RealTypes.h"

namespace EvoEMD {

class Base_Parameter {
protected:
    std::string name;
    REAL value;
    bool updated;
    std::set<Base_Parameter *> base_parameters;  // Parameter in this set depends on current parameter

public:
    Base_Parameter(std::string par_name) : name(par_name), value(0), updated(true){};
    virtual ~Base_Parameter(){};

    bool Is_Independent() { return (base_parameters.size() == 0); }
    void Claim_Dependence(Base_Parameter *par) { base_parameters.insert(par); }
    void Notify() {
        updated = false;
        for (auto &&bp : base_parameters) {
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

class Free_Parameter : public Base_Parameter {
public:
    Free_Parameter(std::string name, REAL default_value = 0) : Base_Parameter(name) {
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
    typedef std::map<std::string, Base_Parameter *> Parameter_List;
    Parameter_List PL;

public:
    static Parameter_Factory &Get_Parameter_Factory();

    void Register_Parameter(Base_Parameter *par);
    bool Set_Parameter(std::string name, REAL value);
    REAL Get_Parameter_Value(std::string name, REAL default_value = 0);
};
}  // namespace EvoEMD
#endif  //_PHYSICS_PARAMETERS_H_
