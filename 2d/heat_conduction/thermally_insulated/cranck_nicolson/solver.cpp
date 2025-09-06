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

using namespace HeatConduction::ThermallyInsulated::CranckNicolson;

Solver::Solver(
    double alpha, double l, double h, 
    double dx, double dy, double r, double epsilon,
    double endTime, double outputTimeStep, string dir
) {
    setALPHA(alpha);

    setL(l);
    setH(h);

    setDX(dx);
    setDY(dy);

    setR(r);

    setEpsilon(epsilon);

    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);

    setDir(dir);
}

double Solver::getALPHA() {
    return alpha;
}

void Solver::setALPHA(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("ALPHA should be > 0, but it is {}", val)
        );
    }

    alpha = val;
}

double Solver::getL() {
    return l;
}

void Solver::setL(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("L should be > 0, but it is {}", val)
        );
    }

    l = val;
}

double Solver::getH() {
    return h;
}

void Solver::setH(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("H should be > 0, but it is {}", val)
        );
    }

    h = val;
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

double Solver::getR() {
    return r;
}

void Solver::setR(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("r should be > 0, but it is {}", val)
        );
    }

    r = val;
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

double Solver::getEndTime() {
    return endTime;
}

void Solver::setEndTime(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("endTime should be > 0, but it is {}", val)
        );
    }

    endTime = val;
}

double Solver::getOutputTimeStep() {
    return outputTimeStep;
}

void Solver::setOutputTimeStep(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("outputTimeStep should be > 0, but it is {}", val)
        );
    }

    outputTimeStep = val;
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
        format("{}/alpha={}, l={}, h={}/dx={}, dy={}, r={}, epsilon={}", dir, alpha, l, h, dx, dy, r, epsilon)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(int nx, int ny, vector<double> &T) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = i*dx;
            double y = j*dy;

            int k = i + j*nx;

            T[k] = cos(M_PI*x/l) * cos(M_PI*y/h);
        }
    }
}

tuple<
    vector<int>, 
    vector<int>, 
    vector<double>
> Solver::buildMatrix(
    int nx, int ny, 
    double rx, double ry
) {
    vector<int> ip(nx*ny + 1);
    vector<int> jp;
    vector<double> e;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int k = i + j*nx;

            ip[k+1] = ip[k] + 1;

            jp.push_back(k);

            if ((i == 0 || i == nx-1) && (j == 0 || j == ny-1)) {
                e.push_back((1 + rx + ry)/4);
            } else if (i == 0 || i == nx-1 || j == 0 || j == ny-1) {
                e.push_back((1 + rx + ry)/2);
            } else {
                e.push_back(1 + rx + ry);
            }            

            if (i > 0) {
                int kl = i - 1 + j*nx;

                ip[k+1] += 1;

                jp.push_back(kl);

                if (j == 0 || j == ny-1) {
                    e.push_back(-rx/4);
                } else {
                    e.push_back(-rx/2);
                }                
            }

            if (j > 0) {
                int kd = i + (j - 1)*nx;

                ip[k+1] += 1;

                jp.push_back(kd);

                if (i == 0 || i == nx-1) {
                    e.push_back(-ry/4);
                } else {
                    e.push_back(-ry/2);
                }
            }
        }
    }

    return {ip, jp, e};
}

void Solver::calculateResidualElements(
    int nx, int ny, 
    double rx, double ry, 
    vector<double> &T, 
    vector<double> &d
) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int k = i + j*nx;

            double Tl = 0;
            double Tr = 0;

            if (i == 0) {
                Tr = T[i + 1 + j*nx];
                Tl = Tr;
            } else if (i == nx-1) {
                Tl = T[i - 1 + j*nx];
                Tr = Tl;
            } else {
                Tl = T[i - 1 + j*nx];
                Tr = T[i + 1 + j*nx];
            }

            double Td = 0;
            double Tu = 0;

            if (j == 0) {
                Tu = T[i + (j + 1)*nx];
                Td = Tu;
            } else if (j == ny-1) {
                Td = T[i + (j - 1)*nx];
                Tu = Td;
            } else {
                Td = T[i + (j - 1)*nx];
                Tu = T[i + (j + 1)*nx];
            }

            d[k] = (1 - rx - ry)*T[k] + rx*(Tr + Tl)/2 + ry*(Tu + Td)/2;     
            
            if ((i == 0 || i == nx-1) && (j == 0 || j == ny-1)) {
                d[k] /= 4;
            } else if (i == 0 || i == nx-1 || j == 0 || j == ny-1) {
                d[k] /= 2;
            }            
        }
    }
}

void Solver::writeData(vector<double> &T, double t, int nx, int ny, int m, path outDir) {
    auto filePath = path(
        format("data.{:03}.vtk", m)
    );

    filePath = outDir / filePath;

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", t) << endl;
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

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx*ny) << endl;
    file << "SCALARS T float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int k = i + j*nx;

            file << format("{:.3f}", T[k]) << endl;
        }
    }
}

void Solver::writeStatistics(vector<tuple<int, double, double>> &statistics, path outDir) {
    auto dirPath = outDir / path("statistics");

    create_directory(dirPath);

    auto filePath = dirPath / path("convergence.csv");

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "n,  t,   max_diff" << endl;

    for (auto t: statistics) {
        file << format(
            "{},  {:.3f},  {:.5f}", 
            get<0>(t),
            get<1>(t),
            get<2>(t)
        ) << endl;
    }
}

void Solver::solve() {
    auto outDir = createDirectory();

    double dt = min(r*dx*dx/alpha, r*dy*dy/alpha);

    double rx = alpha*dt/dx/dx;
    double ry = alpha*dt/dy/dy;

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(h/dy)
    ) + 1;

    vector<double> T(nx*ny);
    vector<double> Tp(nx*ny);

    setInitialCondition(nx, ny, T);

    update(1, T, Tp);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(T, t, nx, ny, m, outDir);

    m += 1;

    vector<int> ip;
    vector<int> jp;
    vector<double> e;
    vector<double> d(nx*ny);

    tie(ip, jp, e) = buildMatrix(nx, ny, rx, ry);

    while (t <= endTime) {         
        calculateResidualElements(nx, ny, rx, ry, T, d);

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
            [&statistics, n](int iter, double error) {
                cout << format("{}: Iteration {}, error {}", n, iter, error) << endl;
            }, 
            ip, jp, e, d, T, epsilon
        );

        t += dt;

        if (t >= tn) {
            writeData(T, t, nx, ny, m, outDir);

            axpby(1, -1, T, Tp, Tp);

            auto maxDiff = maxAbsElem(Tp);

            cout << format("Write data in file t={:.3f}, convergence of T={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            update(1, T, Tp);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}
