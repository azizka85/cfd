#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace PoiseuilleFlow::FreeSurface::PeacemanRachford;

Solver::Solver(
    double nu, double dx, double dy, double r, 
    double endTime, double outputTimeStep, string dir
) {
    setNU(nu);

    setDX(dx);
    setDY(dy);

    setR(r);

    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);

    setDir(dir);
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
        format("{}/nu={}/dx={}, dy={}, r={}", dir, nu, dx, dy, r)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(
    int nx, int ny, 
    vector<vector<double>> &w
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = i*dx;
            double y = j*dy;

            w[i][j] = sin(M_PI*x) * sin(M_PI*y/2);
        }
    }
}

void Solver::setBoundaryCondition(
    int nx, int ny, 
    vector<vector<double>> &w
) {
    for (int i = 0; i < nx; i++) {
        w[i][0] = 0;
    }

    for (int j = 0; j < ny; j++) {
        w[0][j] = 0;
        w[nx-1][j] = 0;
    }
}

double Solver::maxAbsDifference(
    int nx, int ny, 
    vector<vector<double>> &w, 
    vector<vector<double>> &w1
) {
    double maxDiff = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            maxDiff = max(
                maxDiff,
                abs(w[i][j] - w1[i][j])
            );
        }
    }

    return maxDiff;
}

void Solver::updateData(
    int nx, int ny, 
    vector<vector<double>> &w, 
    vector<vector<double>> &w1
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            w[i][j] = w1[i][j];
        }
    }
}

void Solver::writeData(
    vector<vector<double>> &w, 
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
    file << "SCALARS w float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f}", w[i][j]) << endl;
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

    double dt = min(r*dx*dx/nu, r*dy*dy/nu);

    double rx = nu*dt/dx/dx;
    double ry = nu*dt/dy/dy;

    int nx = static_cast<int>(
        ceil(1/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(1/dy)
    ) + 1;

    vector<vector<double>> w(nx, vector<double>(ny));
    vector<vector<double>> w1(nx, vector<double>(ny));
    vector<vector<double>> wp(nx, vector<double>(ny));

    setInitialCondition(nx, ny, w);
    setBoundaryCondition(nx, ny, w);

    updateData(nx, ny, wp, w);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(w, t, nx, ny, m, outDir);

    m += 1;

    vector<double> vx(nx);
    vector<double> vy(ny);
    vector<double> wx(nx);

    while (t <= endTime) {
        setBoundaryCondition(nx, ny, w1); 

        for (int j = 1; j < ny; j++) {
            vx[0] = w[0][j];

            if (j < ny-1) {
                for (int i = 1; i < nx-1; i++) {
                    vx[i] = 0.5*ry*w[i][j+1] + (1 - ry)*w[i][j] + 0.5*ry*w[i][j-1];
                }
            } else {
                for (int i = 1; i < nx-1; i++) {
                    vx[i] = (1 - ry)*w[i][j] + ry*w[i][j-1];
                }
            }

            vx[nx-1] = w[nx-1][j];

            Tridiagonal::solve(-0.5*rx, 0, 1 + rx, 1, 1, -0.5*rx, 0, vx, wx);

            for (int i = 0; i < nx; i++) {
                w1[i][j] = wx[i];
            }
        }

        for (int i = 1; i < nx-1; i++) {
            vy[0] = w1[i][0];

            for (int j = 1; j < ny-1; j++) {
                vy[j] = 0.5*rx*w1[i+1][j] + (1 - rx)*w1[i][j] + 0.5*rx*w1[i-1][j];
            }

            vy[ny-1] = w1[i][ny-1];

            Tridiagonal::solve(-0.5*ry, -ry, 1 + ry, 1, 1 + ry, -0.5*ry, 0, vy, w[i]);
        }

        setBoundaryCondition(nx, ny, w); 

        t += dt;

        if (t >= tn) {
            writeData(w, t, nx, ny, m, outDir);

            auto maxDiff = maxAbsDifference(nx, ny, w, wp);

            cout << format("Write data in file t={:.3f}, convergence of w={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            updateData(nx, ny, wp, w);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}
