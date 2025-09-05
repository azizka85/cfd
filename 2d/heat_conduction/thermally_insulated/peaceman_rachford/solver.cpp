#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace HeatConduction::ThermallyInsulated::PeacemanRachford;

Solver::Solver(
    double alpha, double l, double h,
    double dx, double dy, double r, 
    double endTime, double outputTimeStep, string dir
) {
    setALPHA(alpha);

    setL(l);
    setH(h);

    setDX(dx);
    setDY(dy);

    setR(r);

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
        format("{}/alpha={}, l={}, h={}/dx={}, dy={}, r={}", dir, alpha, l, h, dx, dy, r)
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
            double x = i*dx;
            double y = j*dy;

            T[i][j] = cos(M_PI*x/l) * cos(M_PI*y/h);
        }
    }
}

double Solver::maxAbsDifference(
    int nx, int ny, 
    vector<vector<double>> &T, 
    vector<vector<double>> &T1
) {
    double maxDiff = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            maxDiff = max(
                maxDiff,
                abs(T[i][j] - T1[i][j])
            );
        }
    }

    return maxDiff;
}

void Solver::updateData(
    int nx, int ny, 
    vector<vector<double>> &T, 
    vector<vector<double>> &T1
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            T[i][j] = T1[i][j];
        }
    }
}

void Solver::writeData(
    vector<vector<double>> &T, 
    double t, int nx, int ny, 
    int m, path outDir
) {
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
            file << format("{:.3f}", T[i][j]) << endl;
        }
    }
}

void Solver::writeStatistics(
    vector<tuple<int, double, double>> &statistics, 
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

    vector<vector<double>> T(nx, vector<double>(ny));
    vector<vector<double>> T1(nx, vector<double>(ny));
    vector<vector<double>> Tp(nx, vector<double>(ny));

    setInitialCondition(nx, ny, T);    
    updateData(nx, ny, Tp, T);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(T, t, nx, ny, m, outDir);

    m += 1;

    vector<double> vx(nx);
    vector<double> vy(ny);
    vector<double> Tx(nx);

    while (t <= endTime) {
        for (int i = 0; i < nx; i++) {
            vx[i] = ry*T[i][1] + (1 - ry)*T[i][0];
        }

        Tridiagonal::solve(-0.5*rx, -rx, 1 + rx, 1 + rx, 1 + rx, -0.5*rx, -rx, vx, Tx);

        for (int i = 0; i < nx; i++) {
            T1[i][0] = Tx[i];
        }

        for (int j = 1; j < ny-1; j++) {
            for (int i = 0; i < nx; i++) {
                vx[i] = 0.5*ry*T[i][j+1] + (1 - ry)*T[i][j] + 0.5*ry*T[i][j-1];
            }

            Tridiagonal::solve(-0.5*rx, -rx, 1 + rx, 1 + rx, 1 + rx, -0.5*rx, -rx, vx, Tx);

            for (int i = 0; i < nx; i++) {
                T1[i][j] = Tx[i];
            }
        }

        for (int i = 0; i < nx; i++) {
            vx[i] = (1 - ry)*T[i][ny-1] + ry*T[i][ny-2];
        }

        Tridiagonal::solve(-0.5*rx, -rx, 1 + rx, 1 + rx, 1 + rx, -0.5*rx, -rx, vx, Tx);

        for (int i = 0; i < nx; i++) {
            T1[i][ny-1] = Tx[i];
        }

        for (int j = 0; j < ny; j++) {
            vy[j] = rx*T1[1][j] + (1 - rx)*T1[0][j];
        }

        Tridiagonal::solve(-0.5*ry, -ry, 1 + ry, 1 + ry, 1 + ry, -0.5*ry, -ry, vy, T[0]);

        for (int i = 1; i < nx-1; i++) {
            for (int j = 0; j < ny; j++) {
                vy[j] = 0.5*rx*T1[i+1][j] + (1 - rx)*T1[i][j] + 0.5*rx*T1[i-1][j];
            }

            Tridiagonal::solve(-0.5*ry, -ry, 1 + ry, 1 + ry, 1 + ry, -0.5*ry, -ry, vy, T[i]);
        }

        for (int j = 0; j < ny; j++) {
            vy[j] = (1 - rx)*T1[nx-1][j] + rx*T1[nx-2][j];
        }

        Tridiagonal::solve(-0.5*ry, -ry, 1 + ry, 1 + ry, 1 + ry, -0.5*ry, -ry, vy, T[nx-1]);

        t += dt;

        if (t >= tn) {
            writeData(T, t, nx, ny, m, outDir);

            auto maxDiff = maxAbsDifference(nx, ny, T, Tp);

            cout << format("Write data in file t={:.3f}, convergence of w={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            updateData(nx, ny, Tp, T);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}
