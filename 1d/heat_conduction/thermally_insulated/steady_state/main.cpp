#include <iostream>

#include "solver.h"

using namespace HeatConduction::ThermallyInsulated::SteadyState;

int main() {
    double a = 1;
    double L = 1;
    double dx = 0.01;

    string dir = "data";

    try {
        Solver solver(a, L, dx, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
