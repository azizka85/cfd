#ifndef POISEUILLE_FLOW_CRANCK_NICOLSON_SOLVER_H
#define POISEUILLE_FLOW_CRANCK_NICOLSON_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace PoiseuilleFlow::CranckNicolson {
    class Solver {
        private:
            double nu;
            double dx;
            double dy;
            double r;
            double epsilon;
            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            void setInitialCondition(int nx, int ny,  vector<double>& w);
            void setBoundaryCondition(int nx, int ny,  vector<double>& w);                  

            tuple<vector<int>, vector<int>, vector<double>> buildMatrix(int nx, int ny, double rx, double ry);
            void calculateResidualElements(int nx, int ny, double rx, double ry, vector<double>& w, vector<double>& d);

            void writeData(vector<double>& w, double t, int nx, int ny, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

        public:
            Solver(double nu, double dx, double dy, double r, double epsilon, double endTime, double outputTimeStep, string dir);

            double getNU();
            void setNU(double val);

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
