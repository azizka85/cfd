#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace CouetteFlow::FreeSurface::CranckNicolsonWithTimeVariation::FirstOrderBC;

Solver::Solver(double uBottom, double nu, double h, double dy, double r, double b, double endTime, double outputTimeStep, string dir) {
    setUBottom(uBottom);
    setNU(nu);
    setH(h);
    setDY(dy);
    setR(r);
    setB(b);
    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);
    setDir(dir);
}

double Solver::getUBottom() {
    return uBottom;
}

void Solver::setUBottom(double val) {
    uBottom = val;
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
        format("{}/nu={}, U={}, H={}/dy={}, r={}, b={}", dir, nu, uBottom, h, dy, r, b)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(int ny, vector<double>& u) {
    for (int i = 0; i < ny; i++) {
        u[i] = 0.;
    }
}

void Solver::setBoundaryCondition(int ny, vector<double> &u) {
    if (ny > 0) {
        u[0] = uBottom;
    }
}

void Solver::calculateResidualElements(int ny, double r, vector<double> &u, vector<double> &d) {
    if (ny > 0) {
        d[0] = u[0];
        d[ny-1] = 0;

        for (int i = 1; i < ny - 1; i++) {
            d[i] = r*u[i+1]/2 + (1 - r)*u[i] + r*u[i-1]/2;
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

double Solver::maxAbsDifference(int ny, vector<double> &u, vector<double> &u1) {
    double maxDiff = 0.;

    for (int i = 0; i < ny; i++) {
        maxDiff = max(
            maxDiff,
            abs(u[i] - u1[i])
        );
    }

    return maxDiff;
}

void Solver::updateData(int ny, vector<double> &u, vector<double> &u1) {
    for (int i = 0; i < ny; i++) {
        u[i] = u1[i];
    }
}

void Solver::writeData(vector<double> &u, double t, double dy, int ny, int m, path outDir) {
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
    file << format("DIMENSIONS 1 {} 1", ny) << endl;    
    file << format("POINTS {} float", ny) << endl;

    for (int i = 0; i < ny; i++) {
        double y = dy * i;

        file << format("0.0 {:.3f} 0.0", y) << endl;
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", ny) << endl;
    file << "SCALARS u float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < ny; i++) {
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

    double cr = r;
    double dt = cr*dy*dy/nu;

    int ny = static_cast<int>(
        ceil(h/dy)
    ) + 1;

    vector<double> u(ny);    
    vector<double> up(ny);

    setInitialCondition(ny, u);
    setBoundaryCondition(ny, u);

    updateData(ny, up, u);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(u, t, dy, ny, m, outDir);

    m += 1;

    vector<double> d(ny);

    while (t <= endTime) {
        double al = -cr/2;
        double al1 = -1;

        double ac = 1 + cr;
        double ac0 = 1;
        double ac1 = 1;

        double ar = -cr/2;
        double ar0 = 0;

        calculateResidualElements(ny, cr, u, d);

        Tridiagonal::solve(al, al1, ac, ac0, ac1, ar, ar0, d, u);

        t += dt;

        double cdt = adjustTimeStep(t, dt);

        cr = cr * cdt / dt;
        dt = cdt;

        if (t >= tn) {
            writeData(u, t, dy, ny, m, outDir);

            auto maxDiff = maxAbsDifference(ny, u, up);

            cout << format("Write data in file t={:.3f}, convergence of u={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            updateData(ny, up, u);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);

    cout << "Total number of iterations: " << n;
}
