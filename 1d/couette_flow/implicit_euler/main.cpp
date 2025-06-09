#include <iostream>

#include "solver.h"

using namespace std;

using namespace CouetteFlow::ImplicitEuler;

int main() {
    double uTop = 1;
    double nu = 1;
    double h = 1;
    double dy = 0.01;
    double r = 0.3;
    double endTime = 0.3;
    double outputTimeStep = 0.003;

    string dir = "data";

    try {
        Solver solver(uTop, nu, h, dy, r, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
