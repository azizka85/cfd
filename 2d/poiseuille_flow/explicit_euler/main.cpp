#include <iostream>

#include "solver.h"

using namespace std;

using namespace PoiseuilleFlow::ExplicitEuler;

int main() {
    double nu = 1;
    double dx = 0.01;
    double dy = 0.01;
    double r = 0.15;
    double endTime = 0.3;
    double outputTimeStep = 0.003;

    string dir = "data";

    try {
        Solver solver(nu, dx, dy, r, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
