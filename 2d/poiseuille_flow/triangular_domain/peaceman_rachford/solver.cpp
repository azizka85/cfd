#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace PoiseuilleFlow::TriangularDomain::PeacemanRachford;

Solver::Solver(
    double nu, double a, double dx, double dy, double r, 
    double endTime, double outputTimeStep, string dir
) {
    setNU(nu);

    setA(a);

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

tuple<double, double, double, double> Solver::generateDomain() {
    return {-a/2, -a, a, a};
}

tuple<bool, bool> Solver::pointInDomain(double x, double y) {
    if (isPointInDomain(x, y)) {
        double xl = x - dx;
        double xr = x + dx;

        double yd = y - dy;
        double yu = y + dy;

        if(
            !isPointInDomain(xl, y) ||
            !isPointInDomain(xr, y) ||
            !isPointInDomain(x, yd) ||
            !isPointInDomain(x, yu)
        ) {
            return {true, true};
        }

        return {true, false};
    }

    return {false, false};    
}

bool Solver::isPointInDomain(double x, double y) {
    return y <= 0 && y >= -a*sqrt(3)/2 + sqrt(3)*x && y >= -a*sqrt(3)/2 - sqrt(3)*x;
}

path Solver::createDirectory() {
    auto dirPath = path(
        format("{}/nu={}, a={}/dx={}, dy={}, r={}", dir, nu, a, dx, dy, r)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(
    int nx, int ny, 
    double x0, double y0,
    vector<vector<double>> &w
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = x0 + i*dx;
            double y = y0 + j*dy;

            auto [isDomain, _] = pointInDomain(x, y);

            if (isDomain) {
                w[i][j] = cos(M_PI*x) * sin(M_PI*y);
            } else {
                w[i][j] = 0;
            }
        }
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
    double x0, double y0,
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
            double x = x0 + dx * i;
            double y = y0 + dy * j;

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

    auto [x0, y0, l, h] = generateDomain();

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(h/dy)
    ) + 1;

    vector<vector<double>> w(nx, vector<double>(ny));
    vector<vector<double>> w1(nx, vector<double>(ny));
    vector<vector<double>> wp(nx, vector<double>(ny));

    setInitialCondition(nx, ny, x0, y0, w);

    updateData(nx, ny, wp, w);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(w, t, nx, ny, x0, y0, m, outDir);

    m += 1;        

    while (t <= endTime) {
        for (int i = 0; i < nx; i++) {
            auto [isDomain, _] = pointInDomain(x0 + i*dx, y0 + (ny - 1)*dy);

            if (isDomain) {
                double x = x0 + i*dx;
                double y = y0 + (ny - 1)*dy;

                double wb = cos(M_PI*x)*sin(M_PI*y)*exp(-2*M_PI*M_PI*nu*t);
                double w1b = cos(M_PI*x)*sin(M_PI*y)*exp(-2*M_PI*M_PI*nu*(t + dt));

                w1[i][ny-1] = (wb + w1b)/2;
            }
        }

        for (int j = 0; j < ny-1; j++) {
            vector<double> vx;

            for (int i = 0; i < nx; i++) {
                auto [isDomain, isBoundary] = pointInDomain(x0 + i*dx, y0 + j*dy);

                if (isDomain) {
                    if (isBoundary) {
                        double x = x0 + i*dx;
                        double y = y0 + j*dy;

                        double wb = cos(M_PI*x)*sin(M_PI*y)*exp(-2*M_PI*M_PI*nu*t);
                        double w1b = cos(M_PI*x)*sin(M_PI*y)*exp(-2*M_PI*M_PI*nu*(t + dt));

                        vx.push_back((wb + w1b)/2);
                    } else {
                        vx.push_back(0.5*ry*w[i][j+1] + (1 - ry)*w[i][j] + 0.5*ry*w[i][j-1]);
                    }
                }
            }

            vector<double> wx(vx.size());

            Tridiagonal::solve(-0.5*rx, 0, 1 + rx, 1, 1, -0.5*rx, 0, vx, wx);

            int k = 0;

            for (int i = 0; i < nx; i++) {
                auto [isDomain, _] = pointInDomain(x0 + i*dx, y0 + j*dy);

                if (isDomain) {
                    w1[i][j] = wx[k];

                    k += 1;
                }
            }
        }

        for (int i = 0; i < nx; i++) {
            vector<double> vy;

            for (int j = 0; j < ny; j++) {
                auto [isDomain, isBoundary] = pointInDomain(x0 + i*dx, y0 + j*dy);

                if (isDomain) {
                    if (isBoundary) {
                        double x = x0 + i*dx;
                        double y = y0 + j*dy;

                        vy.push_back(cos(M_PI*x)*sin(M_PI*y)*exp(-2*M_PI*M_PI*nu*(t + dt)));
                    } else {
                        vy.push_back(0.5*rx*w1[i+1][j] + (1 - rx)*w1[i][j] + 0.5*rx*w1[i-1][j]);
                    }
                }
            }

            vector<double> wy(vy.size());

            Tridiagonal::solve(-0.5*ry, 0, 1 + ry, 1, 1, -0.5*ry, 0, vy, wy);

            int k = 0;

            for (int j = 0; j < ny; j++) {
                auto [isDomain, _] = pointInDomain(x0 + i*dx, y0 + j*dy);

                if (isDomain) {
                    w[i][j] = wy[k];

                    k += 1;
                }
            }
        }

        t += dt;

        if (t >= tn) {
            writeData(w, t, nx, ny, x0, y0, m, outDir);

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
