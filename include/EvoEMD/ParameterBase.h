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
    std::set<Parameter_Base *> descendent_parameters;  // Parameter in this set depends on current parameter
    std::set<Parameter_Base *> parent_parameters;

public:
    /**
     * @brief  default ctor, any parameter will be built with a name
     * @note
     * @param  par_name: the name for the parameter
     * @retval
     */
    Parameter_Base(std::string par_name) : name(par_name), value(0), updated(false){};
    virtual ~Parameter_Base(){};

    /**
     * @brief  Checking if current parameter is independent or not
     * @note
     * @retval true for independent, false otherwise
     */
    bool Is_Independent() { return (parent_parameters.size() == 0); }

    /**
     * @brief  Notify that par depends on current parameter
     * @note
     * @param  *par: ptr to the parameter that depends on current parameter
     * @retval None
     */
    void Register_Descendent_Parameter(Parameter_Base *par) { descendent_parameters.insert(par); }

    /**
     * @brief  To notify all parameters that depend on current parameter, that they need to update their values
     * @note
     * @retval None
     */
    void Notify() {
        updated = false;
        for (auto &&bp : descendent_parameters) {
            bp->Notify();
        }
    }

    /**
     * @brief  User use this function to set the value of the parameters
     * @note
     * @param  input: the input value for free parameter, but dummy for derived parameter
     * @retval None
     */
    void Set_Value(REAL input = 0) {
        Notify();
        Update_Value(input);
        updated = true;
    }

    /**
     * @brief  The actual function to update the value of current parameter
     * @note   Any parameter derived from this Base_Parameter should re-implement this function
     * @param  input: the input value for free parameter, but can be ignored for derived parameter
     * @retval None
     */
    virtual void Update_Value(REAL input) = 0;

    /**
     * @brief  Get the value of current parameter
     * @note   This function will take care of updating value if necessary
     * @retval current updated value of the parameter
     */
    REAL Get_Value() {
        if (updated) {
            return value;
        } else {
            Set_Value();
            return value;
        }
    }

    /**
     * @brief  Get the name for the parameter
     * @note
     * @retval
     */
    std::string Get_Name() const { return name; }

    /**
     * @brief  register that current parameter depends on ptr
     * @note
     * @param  *ptr: the ptr to parameter that current parameter depends on
     * @retval None
     */
    void Register_Dependencies(Parameter_Base *ptr) {
        ptr->Register_Descendent_Parameter(this);
        this->parent_parameters.insert(ptr);
    }
    template <typename... Ptrs>
    void Register_Dependencies(Parameter_Base *ptr, Ptrs... params) {
        Register_Dependencies(ptr);
        Register_Dependencies(params...);
    }
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

#define REGISTER_PARAMETER(className) Register_Parameter g_register_parameter_##className(new className)
#define DECLARE_FREE_PARAMETER(paramName, value)                   \
    class param_##paramName : public Free_Parameter {              \
    public:                                                        \
        param_##paramName() : Free_Parameter(#paramName, value){}; \
    };                                                             \
    const Register_Parameter g_register_parameter_##paramName(new param_##paramName)

#define RETRIEVE_PARAMETER(paramName) EvoEMD::Parameter_Factory::Get_Parameter_Factory().Get_Parameter(#paramName)

#define GET_PARAMETER_VALUE(paramName) RETRIEVE_PARAMETER(paramName)->Get_Value()

#endif  //_PARAMETER_BASE_H_
