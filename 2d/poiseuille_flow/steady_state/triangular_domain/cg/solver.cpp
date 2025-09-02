#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/utils.h>
#include <slae/iterative/cg.h>

#include "solver.h"

using namespace SLAE::Utils;
using namespace SLAE::Iterative::CG;

using namespace PoiseuilleFlow::SteadyState::TriangularDomain::CG;

Solver::Solver(
    double g, double nu, double a,
    double dx, double dy, double epsilon,
    string dir
) {
    setG(g);
    setNU(nu);

    setA(a);

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
            format("NU should be > 0, but it is {}", val)
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
            format("A should be > 0, but it is {}", val)
        );
    }

    a = val;
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
    return {-a/2, -a, a, a};
}

tuple<bool, bool> Solver::pointInDomain(double x, double y) {
    if (isPointInDomain(x, y)) {
        double xl = x - dx;
        double xr = x + dx;

        double yd = y - dy;
        double yu = y + dy;

        if(
            !isPointInDomain(xl, y) ||
            !isPointInDomain(xr, y) ||
            !isPointInDomain(x, yd) ||
            !isPointInDomain(x, yu)
        ) {
            return {true, true};
        }

        return {true, false};
    }

    return {false, false};    
}

bool Solver::isPointInDomain(double x, double y) {
    return y <= 0 && y >= -a*sqrt(3)/2 + sqrt(3)*x && y >= -a*sqrt(3)/2 - sqrt(3)*x;
}


path Solver::createDirectory() {
    auto dirPath = path(
        format("{}/g={}, nu={}, a={}/dx={}, dy={}, epsilon={}", dir, g, nu, a, dx, dy, epsilon)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(
    int nx, int ny, 
    vector<double> &w
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {            
            int k = i + j*nx;

            w[k] = 0;
        }
    }
}

tuple<
    vector<int>, 
    vector<int>, 
    vector<double>
> Solver::buildMatrix(
    int nx, int ny, 
    double x0, double y0
) {
    vector<int> ip(nx*ny + 1);
    vector<int> jp;
    vector<double> e;

    double r = dx/dy;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = x0 + i*dx;
            double y = y0 + j*dy;

            int k = i + j*nx;

            auto [isDomain, isBoundary] = pointInDomain(x, y);

            ip[k+1] = ip[k] + 1;

            jp.push_back(k);

            if (!isDomain || isBoundary) {
                e.push_back(1);
            } else {
                e.push_back(-2 - 2*r*r);

                auto [isDomain, isBoundary] = pointInDomain(x - dx, y);

                if (isDomain && !isBoundary) {
                    int kl = i - 1 + j*nx;

                    ip[k+1] += 1;

                    jp.push_back(kl);
                    e.push_back(1);
                }

                tie(isDomain, isBoundary) = pointInDomain(x, y - dy);

                if (isDomain && !isBoundary) {
                    int kd = i + (j - 1)*nx;

                    ip[k+1] += 1;

                    jp.push_back(kd);
                    e.push_back(r*r);
                }
            }                 
        }
    }

    return {ip, jp, e};
}

void Solver::calculateResidualElements(
    int nx, int ny, 
    double x0, double y0, 
    vector<double> &w, 
    vector<double> &d
) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = x0 + i*dx;
            double y = y0 + j*dy;

            auto [isDomain, isBoundary] = pointInDomain(x, y);

            int k = i + j*nx;

            if (isDomain && !isBoundary) {                
                d[k] = dx*dx*g/nu; 
            } else {
                d[k] = 0;
            }           
        }
    }
}

void Solver::writeData(
    vector<double> &w, 
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
            int k = i + j*nx;

            file << format("{:.3f}", w[k]) << endl;
        }
    }
}

void Solver::writeStatistics(vector<tuple<int, double>> &statistics, path outDir) {
    auto dirPath = outDir / path("statistics");

    create_directory(dirPath);

    auto filePath = dirPath / path("convergence.csv");

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "n,  max_diff" << endl;

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

    vector<double> w(nx*ny);
    vector<double> d(nx*ny);

    setInitialCondition(nx, ny, w);        

    vector<int> ip;
    vector<int> jp;
    vector<double> e;

    tie(ip, jp, e) = buildMatrix(nx, ny, x0, y0);    
    calculateResidualElements(nx, ny, x0, y0, w, d);

    vector<tuple<int, double>> statistics;

    statistics.push_back({0, maxAbsElem(d)});

    double tol = epsilon;

    preconditionedCG(
        [&ip, &jp, &e, tol](const vector<double> &r, vector<double> &y) {
            int n = ip.size();

            for (int i = 0; i < n-1; i++) {
                int j = ip[i];
                int k = jp[j];

                double q = abs(e[j]);

                if (q > tol) {
                    y[k] = r[k] / q;
                } else {
                    y[k] = r[k];
                }
            }
        }, 
        [&statistics](int iter, double error) {
            cout << format("Iteration {}, error {}", iter, error) << endl;

            if (error > 1e-7) {
                statistics.push_back({iter, error});
            }
        }, 
        ip, jp, e, d, w, epsilon
    );

    cout << "Write data in file" << endl;

    writeData(w, nx, ny, x0, y0, outDir);
    writeStatistics(statistics, outDir);
}
