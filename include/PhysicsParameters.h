#ifndef _PHYSICS_PARAMETERS_H_
#define _PHYSICS_PARAMETERS_H_

#include <map>
#include <set>
#include <string>

#include "Neutrino.h"

class Base_Parameter {
protected:
    std::string name;
    REAL value;

public:
    Base_Parameter(std::string par_name) : name(par_name){};
    virtual ~Base_Parameter(){};

    virtual void Set_Value(REAL input) { value = input; }
    virtual REAL Get_Value() = 0;
    std::string Get_Name() const { return name; }
};

class Independent_Parameter : public Base_Parameter {
public:
    typedef ROTATION_NUMBER VERSION_TYPE;
    Independent_Parameter(std::string par_name) : Base_Parameter(par_name){};
    ~Independent_Parameter(){};

    virtual void Set_Value(REAL input) {
        value = input;
        ++VERSION_ID;
    }
    virtual REAL Get_Value() { return value; }
    VERSION_TYPE Get_Version_ID() const { return VERSION_ID; }

private:
    VERSION_TYPE VERSION_ID;
};

class Dependent_Parameter : public Base_Parameter {
protected:
    // *
    typedef std::map<std::string, std::pair<Independent_Parameter *, Independent_Parameter::VERSION_TYPE> >
        Independent_Parameter_Set;
    Independent_Parameter_Set IPS;
    bool Is_Updated();
    void Update_Version_ID();
    virtual void Update_Value() = 0;

public:
    Dependent_Parameter(std::string par_name) : Base_Parameter(par_name){};
    ~Dependent_Parameter(){};

    virtual REAL Get_Value();
};

class Parameter_Factory {
private:
    Parameter_Factory();
    ~Parameter_Factory();
    typedef std::map<std::string, Independent_Parameter *> Independent_Parameter_List;
    Independent_Parameter_List IPL;
    typedef std::map<std::string, Dependent_Parameter *> Dependent_Parameter_List;
    Dependent_Parameter_List DPL;

public:
    static Parameter_Factory &Get_Parameter_Factory();

    void Register_Parameter(Independent_Parameter *par);
    void Register_Parameter(Dependent_Parameter *par);
    bool Set_Independent_Parameter(std::string name, REAL value);
    REAL Get_Parameter_Value(std::string name, REAL default_value = 0);
};

#endif  //_PHYSICS_PARAMETERS_H_
