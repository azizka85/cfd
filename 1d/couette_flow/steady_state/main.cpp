#include <iostream>

#include "solver.h"

using namespace CouetteFlow::SteadyState;

int main() {
    double uTop = 1;
    double h = 1;
    double dy = 0.01;

    string dir = "data";

    try {
        Solver solver(uTop, h, dy, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
