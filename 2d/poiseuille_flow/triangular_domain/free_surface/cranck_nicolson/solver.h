#ifndef POISEUILLE_FLOW_TRIANGULAR_DOMAIN_FREE_SURFACE_CRANCK_NICOLSON_SOLVER_H
#define POISEUILLE_FLOW_TRIANGULAR_DOMAIN_FREE_SURFACE_CRANCK_NICOLSON_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace PoiseuilleFlow::TriangularDomain::FreeSurface::CranckNicolson {
    class Solver {
        private:
            double nu;

            double a;

            double dx;
            double dy;

            double r;

            double epsilon;

            double endTime;
            double outputTimeStep;

            string dir;

            tuple<double, double, double, double> generateDomain();
            tuple<bool, bool> pointInDomain(double x, double y);
            bool isPointInDomain(double x, double y);

            path createDirectory();

            void setInitialCondition(int nx, int ny, double x0, double y0,  vector<double>& w);               

            tuple<vector<int>, vector<int>, vector<double>> buildMatrix(int nx, int ny, double x0, double y0, double rx, double ry);
            void calculateResidualElements(int nx, int ny, double x0, double y0, double t, double dt, double rx, double ry, vector<double>& w, vector<double>& d);

            void writeData(vector<double>& w, double t, int nx, int ny, double x0, double y0, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

        public:
            Solver(
                double nu, double a, double dx, double dy, double r, 
                double epsilon, double endTime, double outputTimeStep, string dir
            );

            double getNU();
            void setNU(double val);

            double getA();
            void setA(double val);

            double getDX();
            void setDX(double val);

            double getDY();
            void setDY(double val);

            double getR();
            void setR(double r);

            double getEpsilon();
            void setEpsilon(double val);

            double getEndTime();
            void setEndTime(double val);

            double getOutputTimeStep();
            void setOutputTimeStep(double val);

            string getDir();                                                                        
            void setDir(string val);

            void solve();
    };  
}

#endif
