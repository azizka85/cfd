#include <iostream>

#include "solver.h"

using namespace std;

using namespace CouetteFlow::FreeSurface::CranckNicolsonWithTimeVariation::SecondOrderBC;

int main() {
    double uBottom = 1;
    double nu = 1;
    double h = 1;
    double dy = 0.01;
    double r = 0.5;
    double b = 1.1;
    double endTime = 0.3;
    double outputTimeStep = 0.003;

    string dir = "data";

    try {
        Solver solver(uBottom, nu, h, dy, r, b, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
