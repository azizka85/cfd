#include <format>

#include <fstream>
#include <iostream>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace CouetteFlow::SteadyState;

Solver::Solver(double uTop, double h, double dy, string dir) {
    setUTop(uTop);    
    setH(h);
    setDY(dy);
    setDir(dir);
}

double Solver::getUTop() {
    return uTop;
}

void Solver::setUTop(double val) {
    uTop = val;
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
        format("{}/U={}, H={}/dy={}", dir, uTop, h, dy)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setBoundaryCondition(int ny, vector<double> &u) {
    if (ny > 0) {
        u[0] = 0.;
        u[ny - 1] = uTop;
    }
}

void Solver::writeData(vector<double> &u, double dy, int ny, path outDir) {
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
    file << format("DIMENSIONS 1 {} 1", ny) << endl;    
    file << format("POINTS {} float", ny) << endl;

    for (int i = 0; i < ny; i++) {
        double y = dy * i;

        file << format("0.0 {:.3f} 0.0", y) << endl;
    }
    
    file << format("POINT_DATA {}", ny) << endl;
    file << "SCALARS u float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < ny; i++) {
        file << format("{:.3f}", u[i]) << endl;
    }
}

void Solver::solve() {
    auto outDir = createDirectory();

    int ny = static_cast<int>(
        ceil(h/dy)
    ) + 1;

    vector<double> u(ny);    

    setBoundaryCondition(ny, u);

    double l = 1;
    double l1 = 0;

    double c = -2;
    double c0 = 1;
    double c1 = 1;

    double r = 1;
    double r0 = 0;

    vector<double> u1(ny);

    Tridiagonal::solve(l, l1, c, c0, c1, r, r0, u, u1);

    writeData(u1, dy, ny, outDir);
}
