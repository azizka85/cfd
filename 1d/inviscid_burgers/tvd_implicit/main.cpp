#include <iostream>

#include "solver.h"

using namespace std;

using namespace InviscidBurgers::TVD::Implicit;

int main() {
    double l = 1;
    double dx = 0.01;
    double a = 3;
    double epsilon = 0.1;
    double endTime = 0.3;
    double outputTimeStep = 0.003;

    string dir = "data";

    try {
        Solver solver(l, dx, a, epsilon, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
