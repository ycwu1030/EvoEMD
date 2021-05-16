#include "Neutrino.h"
#include "spdlog_wrapper.h"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    spdlog::set_pattern("[%Y/%m/%d %H:%M:%S][%l][%s:%#:%!] %v");
    Nu_TypeI_SeeSaw nu_model;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << "  " << nu_model.Get_UPMNSij(i, j);
        }
        cout << endl;
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << "  " << nu_model.Get_Yij(i, j);
        }
        cout << endl;
    }

    return 0;
}
