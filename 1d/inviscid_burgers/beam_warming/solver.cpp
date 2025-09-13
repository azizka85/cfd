#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/utils.h>
#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Utils;
using namespace SLAE::Direct;
using namespace InviscidBurgers::BeamWarming;

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

void Solver::calculateRHS(int nx, double dt, vector<double>& u, vector<double>& d) {
    for (int i = 0; i < nx; i++) {
        d[i] = u[i];
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
    vector<double> v(nx);
    vector<double> d(nx);
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

    vector<double> bl(nx);
    vector<double> br(nx);
    vector<double> bc(nx);

    while (t <= endTime) {
        bl[0] = -dt*u[nx-1]/4/dx;
        bl[nx-1] = -dt*u[nx-2]/4/dx;

        br[0] = dt*u[1]/4/dx;
        br[nx-1] = dt*u[0]/4/dx;

        bc[0] = 1;
        bc[nx-1] = 1;

        v[0] = 1;
        v[nx-1] = br[nx-1];

        for (int i = 1; i < nx-1; i++) {
            bl[i] = -dt*u[i-1]/4/dx;
            br[i] = dt*u[i+1]/4/dx;
            bc[i] = 1;
            v[i] = 0;
        }

        calculateRHS(nx, dt, u, d);

        Tridiagonal::solve(bl, bc, br, v, u);
        Tridiagonal::solve(bl, bc, br, d, v);

        double sv = v[0] + bl[0]*v[nx-1];
        double su = u[0] + bl[0]*u[nx-1];

        double p = sv / (1 + su);

        axpby(1, -p, v, u, u);

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
