#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace HeatConduction::ThermallyInsulated::CranckNicolsonWithTimeVariation;

Solver::Solver(double a0, double a1, double alpha, double L, double dx, double r, double b, double endTime, double outputTimeStep, string dir) {
    setA0(a0);
    setA1(a1);
    setAlpha(alpha);
    setL(L);
    setDX(dx);
    setR(r);
    setB(b);
    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);
    setDir(dir);
}

double Solver::getA0() {
    return a0;
}

void Solver::setA0(double val) {
    a0 = val;
}

double Solver::getA1() {
    return a1;
}

void Solver::setA1(double val) {
    a1 = val;
}

double Solver::getAlpha() {
    return alpha;
}

void Solver::setAlpha(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("Alpha should be > 0, but it is {}", val)
        );
    }

    alpha = val;
}

double Solver::getL() {
    return L;
}

void Solver::setL(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("L should be > 0, but it is {}", val)
        );
    }

    L = val;
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
        format("{}/alpha={}, A0={}, A1={}, L={}/dx={}, r={}, b={}", dir, alpha, a0, a1, L, dx, r, b)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(int nx, vector<double>& T) {
    for (int i = 0; i < nx; i++) {
        auto x = i*dx;

        T[i] = a0 + a1 * cos(M_PI * x / L);
    }
}

void Solver::calculateResidualElements(int nx, double r, vector<double> &T, vector<double> &d) {
    if (nx > 0) {
        d[0] = (1 - r)*T[0] + r*T[1];
        d[nx-1] = (1 - r)*T[nx-1] + r*T[nx-2];

        for (int i = 1; i < nx-1; i++) {
            d[i] = r*T[i+1]/2 + (1 - r)*T[i] + r*T[i-1]/2;
        }
    }
}

double Solver::adjustTimeStep(double t, double dt) {
    if (t >= outputTimeStep) {
        return outputTimeStep;
    }

    if (t < outputTimeStep && t + dt >= outputTimeStep) {
        return outputTimeStep - t;
    }

    return b * dt;
}

double Solver::maxAbsDifference(int nx, vector<double> &T, vector<double> &T1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        maxDiff = max(
            maxDiff,
            abs(T[i] - T1[i])
        );
    }

    return maxDiff;
}

void Solver::updateData(int nx, vector<double> &T, vector<double> &T1) {
    for (int i = 0; i < nx; i++) {
        T[i] = T1[i];
    }
}

void Solver::writeData(vector<double> &T, double t, double dx, int nx, int m, path outDir) {
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
    file << "SCALARS T float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < nx; i++) {
        file << format("{:.3f}", T[i]) << endl;
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

    double cr = r;
    double dt = cr*dx*dx/alpha;

    int nx = static_cast<int>(
        ceil(L/dx)
    ) + 1;

    vector<double> T(nx);    
    vector<double> Tp(nx);

    setInitialCondition(nx, T);

    updateData(nx, Tp, T);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(T, t, dx, nx, m, outDir);

    m += 1;

    vector<double> d(nx);

    while (t <= endTime) {
        double al = -cr/2;
        double al1 = -cr;

        double ac = 1 + cr;
        double ac0 = 1 + cr;
        double ac1 = 1 + cr;

        double ar = -cr/2;
        double ar0 = -cr;

        calculateResidualElements(nx, cr, T, d);

        Tridiagonal::solve(al, al1, ac, ac0, ac1, ar, ar0, d, T);

        t += dt;

        double cdt = adjustTimeStep(t, dt);

        cr = cr * cdt / dt;
        dt = cdt;

        if (t >= tn) {
            writeData(T, t, dx, nx, m, outDir);

            auto maxDiff = maxAbsDifference(nx, T, Tp);

            cout << format("Write data in file t={:.3f}, convergence of T={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            updateData(nx, Tp, T);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);

    cout << "Total number of iterations: " << n;
}
