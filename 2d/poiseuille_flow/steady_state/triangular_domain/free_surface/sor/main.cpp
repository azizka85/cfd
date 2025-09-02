#include <iostream>

#include "solver.h"

using namespace std;

using namespace PoiseuilleFlow::SteadyState::TriangularDomain::FreeSurface::SOR;

int main() {
    double a = 1;

    double dx = 0.01;
    double dy = 0.01;

    double epsilon = 1e-6;

    string dir = "data";

    try {
        Solver solver(a, dx, dy, epsilon, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
