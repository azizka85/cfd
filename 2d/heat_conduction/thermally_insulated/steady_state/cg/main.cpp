#include <iostream>

#include "solver.h"

using namespace std;

using namespace HeatConduction::ThermallyInsulated::SteadyState::CG;

int main() {
    double l = 1.0;
    double h = 1.0;

    double dx = 0.01;
    double dy = 0.01;

    double epsilon = 1e-6;

    string dir = "data";

    try {
        Solver solver(l, h, dx, dy, epsilon, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
