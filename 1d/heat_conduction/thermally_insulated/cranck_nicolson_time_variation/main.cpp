#include <iostream>

#include "solver.h"

using namespace std;

using namespace HeatConduction::ThermallyInsulated::CranckNicolsonWithTimeVariation;

int main() {
    double a0 = 1;
    double a1 = 1;
    double alpha = 1;
    double L = 1;
    double dx = 0.01;
    double r = 0.5;
    double b = 1.1;
    double endTime = 0.3;
    double outputTimeStep = 0.003;

    string dir = "data";

    try {
        Solver solver(a0, a1, alpha, L, dx, r, b, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
