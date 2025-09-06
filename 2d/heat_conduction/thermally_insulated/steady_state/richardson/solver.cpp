#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include "solver.h"

using namespace HeatConduction::ThermallyInsulated::SteadyState::Richardson;

Solver::Solver(
    double l, double h,
    double dx, double dy, double epsilon, 
    string dir
) {
    setL(l);
    setH(h);

    setDX(dx);
    setDY(dy);

    setEpsilon(epsilon);

    setDir(dir);
}

double Solver::getL() {
    return l;
}

void Solver::setL(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("l should be > 0, but it is {}", val)
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
            format("h should be > 0, but it is {}", val)
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
        format("{}/l={}, h={}/dx={}, dy={}, epsilon={}", dir, l, h, dx, dy, epsilon)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(
    int nx, int ny, 
    vector<vector<double>> &T
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            T[i][j] = 0;
        }
    }
}

double Solver::calculateMaxAbsResidualElement(
    int nx, int ny, 
    vector<vector<double>>& T
) {
    double maxRes = 0;

    for (int i = 1; i < nx-1; i++) {
        for (int j = 1; j < ny-1; j++) {
            maxRes = max(
                maxRes,
                abs(
                    (T[i+1][j] - 2*T[i][j] + T[i-1][j])/dx/dx + 
                    (T[i][j+1] - 2*T[i][j] + T[i][j-1])/dy/dy + 
                    M_PI*M_PI*(1/l/l + 1/h/h)*cos(M_PI*i*dx/l)*cos(M_PI*j*dy/h)
                )
            );
        }
    }

    return maxRes;
}

void Solver::writeData(
    vector<vector<double>> &T, 
    int nx, int ny, path outDir
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
            double x = dx * i;
            double y = dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << format("POINT_DATA {}", nx*ny) << endl;
    file << "SCALARS T float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f}", T[i][j]) << endl;
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
        ceil(l/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(h/dy)
    ) + 1;

    vector<vector<double>> T(nx, vector<double>(ny));

    setInitialCondition(nx, ny, T);

    double maxR = calculateMaxAbsResidualElement(nx, ny, T);

    vector<tuple<int, double>> statistics;

    statistics.push_back({0, maxR});

    double b = dx/dy;

    int n = 1;

    while (maxR > epsilon) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double Tl = 0;
                double Tr = 0;

                if (i == 0) {
                    Tr = T[i+1][j];
                    Tl = Tr;
                } else if (i == nx-1) {
                    Tl = T[i-1][j];
                    Tr = Tl;
                } else {
                    Tl = T[i-1][j];
                    Tr = T[i+1][j];
                }

                double Td = 0;
                double Tu = 0;

                if (j == 0) {
                    Tu = T[i][j+1];
                    Td = Tu;
                } else if (j == ny-1) {
                    Td = T[i][j-1];
                    Tu = Td;
                } else {
                    Td = T[i][j-1];
                    Tu = T[i][j+1];
                }

                T[i][j] = (
                    Tr + Tl + b*b*Tu + b*b*Td + 
                    dx*dx*M_PI*M_PI*(1/l/l + 1/h/h)*cos(M_PI*i*dx/l)*cos(M_PI*j*dy/h)
                )/2/(b*b + 1);
            }
        }

        maxR = calculateMaxAbsResidualElement(nx, ny, T);

        cout << format(
            "Iteration={}, convergence of w={:.6f}", 
            n, maxR
        ) << endl;

        statistics.push_back({n, maxR});

        n += 1;
    }

    writeData(T, nx, ny, outDir);
    writeStatistics(statistics, outDir);
}
