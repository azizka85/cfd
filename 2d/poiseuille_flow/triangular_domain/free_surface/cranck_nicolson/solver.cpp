#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/utils.h>
#include <slae/iterative/cg.h>

#include "solver.h"

using namespace SLAE::Utils;
using namespace SLAE::Iterative::CG;

using namespace PoiseuilleFlow::TriangularDomain::FreeSurface::CranckNicolson;

Solver::Solver(
    double nu, double a, double dx, double dy, double r, 
    double epsilon, double endTime, double outputTimeStep, string dir
) {
    setNU(nu);

    setA(a);

    setDX(dx);
    setDY(dy);

    setR(r);

    setEpsilon(epsilon);

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
        format("{}/nu={}, a={}/dx={}, dy={}, r={}, epsilon={}", dir, nu, a, dx, dy, r, epsilon)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(
    int nx, int ny, 
    double x0, double y0,
    vector<double> &w
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = x0 + i*dx;
            double y = y0 + j*dy;

            auto [isDomain, _] = pointInDomain(x, y);

            int k = i + j*nx;

            if (isDomain) {
                w[k] = cos(M_PI*x) * cos(M_PI*y);
            } else {
                w[k] = 0;
            }
        }
    }
}

tuple<
    vector<int>, 
    vector<int>, 
    vector<double>
> Solver::buildMatrix(
    int nx, int ny, 
    double x0, double y0,
    double rx, double ry
) {
    vector<int> ip(nx*ny + 1);
    vector<int> jp;
    vector<double> e;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = x0 + i*dx;
            double y = y0 + j*dy;

            int k = i + j*nx;

            auto [isDomain, isBoundary] = pointInDomain(x, y);

            ip[k+1] = ip[k] + 1;

            jp.push_back(k);

            if (!isDomain || (isBoundary && (j < ny-1 || i == 0 || i == nx-1))) {
                e.push_back(1);
            } else {
                if (j < ny-1) {
                    e.push_back(1 + rx + ry);
                } else {
                    e.push_back((1 + rx + ry)/2);
                } 

                auto [isDomain, isBoundary] = pointInDomain(x - dx, y);

                if (isDomain && !isBoundary) {
                    int kl = i - 1 + j*nx;

                    ip[k+1] += 1;

                    jp.push_back(kl);
                    
                    if (j < ny-1) {
                        e.push_back(-rx/2);
                    } else {
                        e.push_back(-rx/4);
                    }
                }

                tie(isDomain, isBoundary) = pointInDomain(x, y - dy);

                if (isDomain && !isBoundary) {
                    int kd = i + (j - 1)*nx;

                    ip[k+1] += 1;

                    jp.push_back(kd);
                    e.push_back(-ry/2);
                }
            }                 
        }
    }

    return {ip, jp, e};
}

void Solver::calculateResidualElements(
    int nx, int ny, 
    double x0, double y0, 
    double t, double dt,
    double rx, double ry, 
    vector<double> &w, 
    vector<double> &d
) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = x0 + i*dx;
            double y = y0 + j*dy;

            auto [isDomain, isBoundary] = pointInDomain(x, y);

            int k = i + j*nx;

            if (isDomain) {                
                if(isBoundary) {
                    if (j == ny-1 && i > 0 && i < nx-1) {
                        int kl = i - 1 + (ny - 1)*nx;
                        int kr = i + 1 + (ny - 1)*nx;
                        int kd = i + (ny - 2)*nx;

                        d[k] = (1 - rx - ry)*w[k]/2 + rx*(w[kr] + w[kl])/4 + ry*w[kd]/2;
                    } else {
                        d[k] = cos(M_PI*x)*cos(M_PI*y)*exp(-2*M_PI*M_PI*nu*t);
                    }                    
                } else {
                    tie(isDomain, isBoundary) = pointInDomain(x - dx, y);

                    double wl = 0;

                    if (isDomain) {
                        if (isBoundary) {
                            wl = cos(M_PI*(x - dx))*cos(M_PI*y)*exp(-2*M_PI*M_PI*nu*t);
                            d[k] += rx*cos(M_PI*(x - dx))*cos(M_PI*y)*exp(-2*M_PI*M_PI*nu*(t + dt))/2;
                        } else {
                            wl = w[i - 1 + j*nx];
                        }
                    }

                    tie(isDomain, isBoundary) = pointInDomain(x + dx, y);

                    double wr = 0;

                    if (isDomain) {
                        if (isBoundary) {
                            wr = cos(M_PI*(x + dx))*cos(M_PI*y)*exp(-2*M_PI*M_PI*nu*t);
                            d[k] += rx*cos(M_PI*(x + dx))*cos(M_PI*y)*exp(-2*M_PI*M_PI*nu*(t + dt))/2;
                        } else {
                            wr = w[i + 1 + j*nx];
                        }
                    }

                    tie(isDomain, isBoundary) = pointInDomain(x, y - dy);

                    double wd = 0;

                    if (isDomain) {
                        if (isBoundary) {
                            wd = cos(M_PI*x)*cos(M_PI*(y - dy))*exp(-2*M_PI*M_PI*nu*t);
                            d[k] += ry*cos(M_PI*x)*cos(M_PI*(y - dy))*exp(-2*M_PI*M_PI*nu*(t + dt))/2;
                        } else {
                            wd = w[i + (j - 1)*nx];
                        }
                    }

                    tie(isDomain, isBoundary) = pointInDomain(x, y + dy);

                    double wu = 0;

                    if (isDomain) {
                        if (isBoundary && (j < ny-2 || i == 0 || i == nx-1)) {
                            wu = cos(M_PI*x)*cos(M_PI*(y + dy))*exp(-2*M_PI*M_PI*nu*t);
                            d[k] += ry*cos(M_PI*x)*cos(M_PI*(y + dy))*exp(-2*M_PI*M_PI*nu*(t + dt))/2;
                        } else {
                            wu = w[i + (j + 1)*nx];
                        }
                    }                           

                    d[k] += (1 - rx - ry)*w[k] + rx*(wr + wl)/2 + ry*(wu + wd)/2;                    
                }
            } else {
                d[k] = 0;
            }           
        }
    }
}

void Solver::writeData(
    vector<double> &w, 
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
            int k = i + j*nx;

            file << format("{:.3f}", w[k]) << endl;
        }
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

    vector<double> w(nx*ny);
    vector<double> wp(nx*ny);

    setInitialCondition(nx, ny, x0, y0, w);    

    update(1, w, wp);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(w, t, nx, ny, x0, y0, m, outDir);

    m += 1;

    vector<int> ip;
    vector<int> jp;
    vector<double> e;
    vector<double> d(nx*ny);

    tie(ip, jp, e) = buildMatrix(nx, ny, x0, y0, rx, ry);

    while (t <= endTime) {         
        calculateResidualElements(nx, ny, x0, y0, t, dt, rx, ry, w, d);

        double tol = epsilon;

        preconditionedCG(
            [&ip, &jp, &e, tol](const vector<double> &r, vector<double> &y) {
                int n = ip.size();

                for (int i = 0; i < n-1; i++) {
                    int j = ip[i];
                    int k = jp[j];

                    double q = abs(e[j]);

                    if (q > tol) {
                        y[k] = r[k] / q;
                    } else {
                        y[k] = r[k];
                    }
                }
            }, 
            [&statistics, n](int iter, double error) {
                cout << format("{}: Iteration {}, error {}", n, iter, error) << endl;
            }, 
            ip, jp, e, d, w, epsilon
        );

        t += dt;

        if (t >= tn) {
            writeData(w, t, nx, ny, x0, y0, m, outDir);

            axpby(1, -1, w, wp, wp);

            auto maxDiff = maxAbsElem(wp);

            cout << format("Write data in file t={:.3f}, convergence of w={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            update(1, w, wp);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}
