#include <format>

#include <fstream>
#include <iostream>

#include <slae/direct/augmented_tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace HeatConduction::ThermallyInsulated::SteadyState;

Solver::Solver(double a, double L, double dx, string dir) {
    setA(a);    
    setL(L);
    setDX(dx);
    setDir(dir);
}

double Solver::getA() {
    return a;
}

void Solver::setA(double val) {
    a = val;
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
        format("{}/A={}, L={}/dx={}", dir, a, L, dx)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::writeData(vector<double> &T, double dx, int nx, path outDir) {
    auto filePath = path("data.vtk");

    filePath = outDir / filePath;

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "# vtk DataFile Version 3.0" << endl;
    file << "TIME 0" << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} 1 1", nx) << endl;    
    file << format("POINTS {} float", nx) << endl;

    for (int i = 0; i < nx; i++) {
        double x = dx * i;

        file << format("{:.3f} 0.0 0.0", x) << endl;
    }
    
    file << format("POINT_DATA {}", nx) << endl;
    file << "SCALARS T float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < nx; i++) {
        file << format("{:.3f}", T[i]) << endl;
    }
}

void Solver::solve() {
    auto outDir = createDirectory();

    int nx = static_cast<int>(
        ceil(L/dx)
    ) + 1;

    vector<double> T(nx);    

    double l = 1;
    
    double c = -2;

    double r = 1;
    double r0 = 2;

    double f = 1;
    double f0 = 0.5;
    
    double h = a*L / dx;

    AugmentedTridiagonal::solve(l, c, c, r, r0, f, f0, f0, h, T);

    writeData(T, dx, nx, outDir);
}
