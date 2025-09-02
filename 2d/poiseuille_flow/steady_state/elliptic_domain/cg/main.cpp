#include <iostream>

#include "solver.h"

using namespace std;

using namespace PoiseuilleFlow::SteadyState::EllipticDomain::CG;

int main() {
    double g = 1;
    double nu = 1;

    double a = 2;
    double b = 1;

    double dx = 0.01;
    double dy = 0.01;

    double epsilon = 1e-6;

    string dir = "data";

    try {
        Solver solver(g, nu, a, b, dx, dy, epsilon, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
