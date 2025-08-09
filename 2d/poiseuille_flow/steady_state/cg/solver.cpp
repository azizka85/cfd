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

using namespace PoiseuilleFlow::SteadyState::CG;

Solver::Solver(double dx, double dy, double epsilon, string dir) {
    setDX(dx);
    setDY(dy);

    setEpsilon(epsilon);

    setDir(dir);
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

path Solver::createDirectory() {
    auto dirPath = path(
        format("{}/dx={}, dy={}, epsilon={}", dir, dx, dy, epsilon)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(int nx, int ny, vector<double> &w) {
    for (int j = 1; j < ny - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            int k = i - 1 + (j - 1)*(nx - 2);

            w[k] = 0;
        }
    }
}

void Solver::setBoundaryCondition(
    int nx, int ny, 
    vector<double> &w
) {
    for (int i = 0; i < nx; i++) {
        int k = i;

        w[k] = 0;

        k = i + (ny - 1)*nx;

        w[k] = 0;
    }

    for (int j = 0; j < ny; j++) {
        int k = j*nx;

        w[k] = 0;

        k = nx - 1 + j*nx;

        w[k] = 0;
    }
}

tuple<vector<int>, vector<int>, vector<double>> Solver::buildMatrix(int nx, int ny) {
    vector<int> ip(nx*ny + 1);
    vector<int> jp;
    vector<double> e;

    double b = dx/dy;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int k = i + j*nx;

            ip[k+1] = ip[k] + 1;

            jp.push_back(k);
            
            if (i == 0 || i == nx-1 || j == 0 || j == ny-1) {
                e.push_back(1);
            } else {
                e.push_back(-2 - 2*b*b);

                if (i > 1) {
                    int kl = i - 1 + j*nx;

                    ip[k+1] += 1;

                    jp.push_back(kl);
                    e.push_back(1);
                }

                if (j > 1) {
                    int kd = i + (j - 1)*nx;

                    ip[k+1] += 1;

                    jp.push_back(kd);
                    e.push_back(b*b);
                }
            }
        }
    }

    return {ip, jp, e};
    
}

void Solver::calculateResidualElements(
    int nx, int ny, 
    vector<double> &w, 
    vector<double> &d
) {
    for (int i = 0; i < nx; i++) {
        int k = i;

        d[k] = w[k];

        k = i + (ny - 1)*nx;

        d[k] = w[k];
    }

    for (int j = 1; j < ny-1; j++) {
        int k = j*nx;

        d[k] = w[k];

        k = nx - 1 + j*nx;

        d[k] = w[k];
    }

    double b = dx/dy;

    for (int j = 1; j < ny-1; j++) {
        for (int i = 1; i < nx-1; i++) {
            int k = i + j*nx;

            d[k] = -2*dx*dx*M_PI*M_PI*sin(M_PI*i*dx)*sin(M_PI*j*dy);

            if (i == 1) {
                int kl = i - 1 + j*nx;

                d[k] -= w[kl];
            }

            if (i == nx-2) {
                int kr = i + 1 + j*nx;

                d[k] -= w[kr];
            }

            if (j == 1) {
                int kd = i + (j - 1)*nx;

                d[k] -= b*b*w[kd];
            }

            if (j == ny-2) {
                int ku = i + (j + 1)*nx;

                d[k] -= b*b*w[ku];
            }
        }
    }
}

void Solver::writeData(
    vector<double> &w, int nx, int ny, path outDir) {
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
            double x = dx * i;
            double y = dy * j;

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

    int nx = static_cast<int>(
        ceil(1/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(1/dy)
    ) + 1;

    vector<double> w(nx*ny);
    vector<double> d(nx*ny);

    setInitialCondition(nx, ny, w); 
    setBoundaryCondition(nx, ny, w);       

    vector<int> ip;
    vector<int> jp;
    vector<double> e;

    tie(ip, jp, e) = buildMatrix(nx, ny);    
    calculateResidualElements(nx, ny, w, d);    

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

            statistics.push_back({iter, error});
        }, 
        ip, jp, e, d, w, epsilon
    );

    cout << "Write data in file" << endl;

    writeData(w, nx, ny, outDir);
    writeStatistics(statistics, outDir);
}
