#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include "solver.h"

using namespace PoiseuilleFlow::SteadyState::EllipticDomain::SOR;

Solver::Solver(
    double g, double nu, 
    double a, double b,
    double dx, double dy, double epsilon, 
    string dir
) {
    setG(g);
    setNU(nu);

    setA(a);
    setB(b);

    setDX(dx);
    setDY(dy);

    setEpsilon(epsilon);

    setDir(dir);
}

double Solver::getG() {
    return g;
}

void Solver::setG(double val) {
    g = val;
}

double Solver::getNU() {
    return nu;
}

void Solver::setNU(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("nu should be > 0, but it is {}", val)
        );
    }

    nu = val;
}

double Solver::getA() {
    return a;
}   

void Solver::setA(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("a should be > 0, but it is {}", val)
        );
    }

    a = val;
}

double Solver::getB() {
    return b;
}

void Solver::setB(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("b should be > 0, but it is {}", val)
        );
    }

    b = val;
}

double Solver::getDX() {
    return dx;
}

void Solver::setDX(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("dx should be > 0, but it is {}", val)
        );
    }

    dx = val;
}

double Solver::getDY() {
    return dy;
}

void Solver::setDY(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("dy should be > 0, but it is {}", val)
        );
    }

    dy = val;
}

double Solver::getEpsilon() {
    return epsilon;
}

void Solver::setEpsilon(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("epsilon should be > 0, but it is {}", val)
        );
    }

    epsilon = val;
}

string Solver::getDir() {
    return dir;
}

void Solver::setDir(string val) {
    if (val.empty()) {
        throw runtime_error(
            format("dir should not be empty, but it is {}", val)
        );
    }

    dir = val;
}

tuple<double, double, double, double> Solver::generateDomain() {
    return {-a, -b, 2*a, 2*b};
}

tuple<bool, bool> Solver::pointInDomain(double x, double y) {
    if (x*x/a/a + y*y/b/b <= 1) {
        double xl = x - dx;
        double xr = x + dx;

        double yd = y - dy;
        double yu = y + dy;

        if (
            xl*xl/a/a + y*y/b/b > 1 ||
            xr*xr/a/a + y*y/b/b > 1 ||
            x*x/a/a + yd*yd/b/b > 1 ||
            x*x/a/a + yu*yu/b/b > 1
        ) {
            return {true, true};
        }

        return {true, false};
    }

    return {false, false};    
}

path Solver::createDirectory() {
    auto dirPath = path(
        format("{}/g={}, nu={}/a={}, b={}/dx={}, dy={}, epsilon={}", dir, g, nu, a, b, dx, dy, epsilon)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(
    int nx, int ny, 
    vector<vector<double>> &w
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            w[i][j] = 0;
        }
    }
}

double Solver::calculateMaxAbsResidualElement(
    int nx, int ny, 
    double x0, double y0,
    vector<vector<double>>& w
) {
    double maxRes = 0;

    for (int i = 1; i < nx-1; i++) {
        for (int j = 1; j < ny-1; j++) {
            double x = x0 + i*dx;
            double y = y0 + j*dy;

            auto [isDomain, isBoundary] = pointInDomain(x, y);

            if (isDomain && !isBoundary) {
                maxRes = max(
                    maxRes,
                    abs(
                        (w[i+1][j] - 2*w[i][j] + w[i-1][j])/dx/dx + 
                        (w[i][j+1] - 2*w[i][j] + w[i][j-1])/dy/dy - 
                        g/nu
                    )
                );
            }
        }
    }

    return maxRes;
}

void Solver::writeData(
    vector<vector<double>> &w, 
    int nx, int ny, 
    double x0, double y0,
    path outDir
) {
    auto filePath = path("data.vtk");

    filePath = outDir / filePath;

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", 0.) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} 1", nx, ny) << endl;    
    file << format("POINTS {} float", nx*ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = x0 + dx * i;
            double y = y0 + dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << format("POINT_DATA {}", nx*ny) << endl;
    file << "SCALARS w float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f}", w[i][j]) << endl;
        }
    }
}

void Solver::writeStatistics(
    vector<tuple<int, double>> &statistics, 
    path outDir
) {
    auto dirPath = outDir / path("statistics");

    create_directory(dirPath);

    auto filePath = dirPath / path("convergence.csv");

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "n,  error" << endl;

    for (auto t: statistics) {
        file << format(
            "{},  {:.7f}", 
            get<0>(t),
            get<1>(t)
        ) << endl;
    }
}

void Solver::solve() {
    auto outDir = createDirectory();

    auto [x0, y0, l, h] = generateDomain();

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(h/dy)
    ) + 1;

    vector<vector<double>> w(nx, vector<double>(ny));

    setInitialCondition(nx, ny, w);

    double maxR = calculateMaxAbsResidualElement(nx, ny, x0, y0, w);

    vector<tuple<int, double>> statistics;

    statistics.push_back({0, maxR});

    double r = dx/dy;
    double f = (cos(M_PI/(nx-1)) + r*r*cos(M_PI/(ny-1)))/(1 + r*r);
    double k = f*f;
    double omega = 2*(1 - sqrt(1 - k))/k;

    int n = 1;

    while (maxR > epsilon) {
        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                auto [isDomain, isBoundary] = pointInDomain(x0 + i*dx, y0 + j*dy);

                if (isDomain && !isBoundary) {
                    w[i][j] += omega*(
                        w[i+1][j] + w[i-1][j] + 
                        r*r*w[i][j+1] + r*r*w[i][j-1] - 
                        2*(1 + r*r)*w[i][j] - 
                        dx*dx*g/nu
                    )/2/(r*r + 1);
                }
            }
        }

        maxR = calculateMaxAbsResidualElement(nx, ny, x0, y0, w);

        cout << format(
            "Iteration={}, convergence of w={:.6f}", 
            n, maxR
        ) << endl;

        statistics.push_back({n, maxR});

        n += 1;
    }

    writeData(w, nx, ny, x0, y0, outDir);
    writeStatistics(statistics, outDir);
}
