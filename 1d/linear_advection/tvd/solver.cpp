#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include "solver.h"

using namespace LinearAdvection::TVD;

Solver::Solver(double u, double d, double l, double dx, double a, double epsilon, double endTime, double outputTimeStep, string dir) {
    setU(u);
    setD(d);
    setL(l);
    setDX(dx);
    setA(a);
    setEpsilon(epsilon);
    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);
    setDir(dir);
}

double Solver::getU() {
    return u;
}

void Solver::setU(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("U should be > 0, but it is {}", val)
        );
    }

    u = val;
}

double Solver::getD() {
    return d;
}

void Solver::setD(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("D should be > 0, but it is {}", val)
        );
    }

    d = val;
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

double Solver::getEpsilon() {
    return epsilon;
}

void Solver::setEpsilon(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("Epsilon should be > 0, but it is {}", val)
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

path Solver::createDirectory() {
    auto dirPath = path(
        format("{}/U={}, d={}, L={}/dx={}, a={}, epsilon={}", dir, u, d, l, dx, a, epsilon)
    );

    create_directories(dirPath);    

    return dirPath;
}

void Solver::setInitialCondition(int nx, vector<double>& C) {
    for (int i = 0; i < nx; i++) {
        double x = dx * i;

        if (x <= d) {
            C[i] = 1;
        } else {
            C[i] = 0;
        }
    }
}

double Solver::maxAbsDifference(int nx, vector<double> &C, vector<double> &C1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        maxDiff = max(
            maxDiff,
            abs(C[i] - C1[i])
        );
    }

    return maxDiff;
}

void Solver::updateData(int nx, vector<double> &C, vector<double> &C1) {
    for (int i = 0; i < nx; i++) {
        C[i] = C1[i];
    }
}

void Solver::writeData(vector<double> &C, double t, double dx, int nx, int m, path outDir) {
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
    file << "SCALARS C float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < nx; i++) {
        file << format("{:.3f}", C[i]) << endl;
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

double Solver::q(double x) {
    if (abs(x) < 2*epsilon) {
        return x*x/4/epsilon + epsilon;
    }

    return abs(x);
}

void Solver::solve() {
    auto outDir = createDirectory();

    double dt = a*dx/u;

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    vector<double> C(nx);
    vector<double> C1(nx);
    vector<double> Cp(nx);

    setInitialCondition(nx, C);

    updateData(nx, Cp, C);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(C, t, dx, nx, m, outDir);

    m += 1;

    while (t <= endTime) {
        for (int i = 0; i < nx; i++) {
            double Cl = C[1];
            double Cr = C[nx-2];

            if (i > 0) {
                Cl = C[i-1];
            }

            if (i < nx-1) {
                Cr = C[i+1];
            }

            double Cll = C1[2];
            double Crr = C1[nx-3];

            if (i == 1) {
                Cll = C1[1];
            } else if (i > 1) {
                Cll = C1[i-2];
            }

            if (i == nx-2) {
                Crr = C1[nx-2];
            } else if (i < nx-2) {
                Crr = C1[i+2];
            }

            double f = u*C[i];
            double fr = u*Cr;
            double fl = u*Cl;

            double fp = (f + fr)/2;
            double fm = (fl + f)/2;

            double ar = a;
            double al = a;

            double arr = a;
            double all = a;

            double dr = Cr - C[i];
            double dl = C[i] - Cl;

            double drr = Crr - Cr;
            double dll = Cl - Cll;

            double gr = (q(ar) - ar*ar)*dr/2;
            double gl = (q(al) - al*al)*dl/2;

            double grr = (q(arr) - arr*arr)*drr/2;
            double gll = (q(all) - all*all)*dll/2;

            double sr = gr > 0 ? 1 : -1;
            double sl = gl > 0 ? 1 : -1;

            double srr = grr > 0 ? 1 : -1;            

            double g = sr*max(0., min(abs(gr), gl*sr));
            double gp = srr*max(0., min(abs(grr), gr*srr));
            double gm = sl*max(0., min(abs(gl), gll*sl));

            double br = dr > 0 ? (gp - g)/dr : 0;
            double bl = dl > 0 ? (g - gm)/dl : 0;

            double hp = (g + gp - q(ar + br)*dr)/2;
            double hm = (g + gm - q(al + bl)*dl)/2;

            C1[i] = C[i] - (fp - fm)*dt/dx - hp + hm;
        }

        updateData(nx, C, C1);   

        t += dt;

        if (t >= tn) {
            writeData(C, t, dx, nx, m, outDir);

            auto maxDiff = maxAbsDifference(nx, C, Cp);

            cout << format("Write data in file t={:.3f}, convergence of u={:.5f}", t, maxDiff) << endl;

            statistics.push_back({n, tn, maxDiff});

            updateData(nx, Cp, C);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}
