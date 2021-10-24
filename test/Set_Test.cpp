#include <iostream>
#include <set>

using namespace std;

int main(int argc, char const *argv[]) {
    set<double *> sd;
    double a = 1;
    double b = 2;
    // double c = 3;
    // double d = 4;

    double *p1 = &a;
    double *p2 = &b;
    double *p3 = &a;
    sd.insert(p1);
    sd.insert(p2);
    sd.insert(p3);

    for (auto &&i : sd) {
        cout << i << "\t" << *i << endl;
    }

    return 0;
}
