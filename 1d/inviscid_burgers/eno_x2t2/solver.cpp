#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include "solver.h"

using namespace InviscidBurgers::ENO::X2T2;

Solver::Solver(double l, double dx, double a, double endTime, double outputTimeStep, string dir) {
    setL(l);
    setDX(dx);
    setA(a);
    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);
    setDir(dir);
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
        format("{}/L={}, dx={}, a={}", dir, l, dx, a)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(int nx, vector<double>& u) {
    for (int i = 0; i < nx; i++) {
        double x = dx * i;

        u[i] = sin(2.*M_PI*x/l);
    }
}

double Solver::maxAbsDifference(int nx, vector<double> &u, vector<double> &u1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        maxDiff = max(
            maxDiff,
            abs(u[i] - u1[i])
        );
    }

    return maxDiff;
}

void Solver::updateData(int nx, vector<double> &u, vector<double> &u1) {
    for (int i = 0; i < nx; i++) {
        u[i] = u1[i];
    }
}

void Solver::writeData(vector<double> &u, double t, double dx, int nx, int m, path outDir) {
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
    file << format("DIMENSIONS {} 1 1", nx) << endl;    
    file << format("POINTS {} float", nx) << endl;

    for (int i = 0; i < nx; i++) {
        double x = dx * i;

        file << format("{:.3f} 0.0 0.0", x) << endl;
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx) << endl;
    file << "SCALARS u float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < nx; i++) {
        file << format("{:.3f}", u[i]) << endl;
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

    double dt = a*dx;

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    vector<double> u(nx);
    vector<double> u1(nx);
    vector<double> up(nx);

    setInitialCondition(nx, u);

    updateData(nx, up, u);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(u, t, dx, nx, m, outDir);

    m += 1;

    while (t <= endTime) {
        for (int i = 0; i < nx; i++) {
            double ul = u[nx-1];
            double ur = u[0];

            if (i > 0) {
                ul = u[i-1];
            }

            if (i < nx-1) {
                ur = u[i+1];
            }

            double ull = u[nx-2];
            double urr = u[1];

            if (i == 1) {
                ull = u[nx-1];
            } else if (i > 1) {
                ull = u[i-2];
            }

            if (i == nx-2) {
                urr = u[0];
            } else if (i < nx-2) {
                urr = u[i+2];
            }

            double vl = (u[i] - ul)/2;
            double vr = (ur - u[i])/2;

            double v = abs(vl) < abs(vr) ? vl : vr;

            double a = ur == u[i] ? u[i] : (ur*ur - u[i]*u[i])/(ur - u[i])/2;

            if (a > 0) {
                double up = u[i] + v;

                vl = (ul - ull)/2;
                vr = (u[i] - ul)/2;

                v = abs(vl) < abs(vr) ? vl : vr;

                double um = ul + v;

                u1[i] = u[i] - dt*(up*up - um*um)/2/dx;
            } else {
                double um = u[i] - v;

                vl = (ur - u[i])/2;
                vr = (urr - ur)/2;

                v = abs(vl) < abs(vr) ? vl : vr;

                double up = ur - v;

                u1[i] = u[i] - dt*(up*up - um*um)/2/dx;
            }
        }

        for (int i = 0; i < nx; i++) {
            double ul = u1[nx-1];
            double ur = u1[0];

            if (i > 0) {
                ul = u1[i-1];
            }

            if (i < nx-1) {
                ur = u1[i+1];
            }

            double ull = u1[nx-2];
            double urr = u1[1];

            if (i == 1) {
                ull = u1[nx-1];
            } else if (i > 1) {
                ull = u1[i-2];
            }

            if (i == nx-2) {
                urr = u1[0];
            } else if (i < nx-2) {
                urr = u1[i+2];
            }

            double vl = (u1[i] - ul)/2;
            double vr = (ur - u1[i])/2;

            double v = abs(vl) < abs(vr) ? vl : vr;

            double a = ur == u1[i] ? u1[i] : (ur*ur - u1[i]*u1[i])/(ur - u1[i])/2;

            if (a > 0) {
                double up = u1[i] + v;

                vl = (ul - ull)/2;
                vr = (u1[i] - ul)/2;

                v = abs(vl) < abs(vr) ? vl : vr;

                double um = ul + v;

                u[i] = (u[i] + u1[i])/2 - dt*(up*up - um*um)/4/dx;
            } else {
                double um = u1[i] - v;

                vl = (ur - u1[i])/2;
                vr = (urr - ur)/2;

                v = abs(vl) < abs(vr) ? vl : vr;

                double up = ur - v;

                u[i] = (u[i] + u1[i])/2 - dt*(up*up - um*um)/4/dx;
            }
        }

        t += dt;

        if (t >= tn) {
            writeData(u, t, dx, nx, m, outDir);

            auto maxDiff = maxAbsDifference(nx, u, up);

            cout << format("Write data in file t={:.3f}, convergence of u={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            updateData(nx, up, u);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}
