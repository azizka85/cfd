#ifndef COUETTE_FLOW_GRAVITY_CRANCK_NICOLSON_TIME_VARIATION_SOLVER_H
#define COUETTE_FLOW_GRAVITY_CRANCK_NICOLSON_TIME_VARIATION_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace CouetteFlow::Gravity::CranckNicolsonWithTimeVariation {
    class Solver {
        private:
            double gx;
            double nu;
            double h;
            double dy;
            double r;
            double b;
            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            void setInitialCondition(int ny, vector<double>& u);
            void setBoundaryCondition(int ny, vector<double>& u);       
            
            void calculateResidualElements(int ny, double dt, double r, vector<double>& u, vector<double>& d);

            double adjustTimeStep(double t, double dt);
            
            double maxAbsDifference(int ny, vector<double>& u, vector<double>& u1);

            void updateData(int ny, vector<double>& u, vector<double>& u1);
            void writeData(vector<double>& u, double t, double dy, int ny, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

        public:
            Solver(double gx, double nu, double h, double dy, double r, double b, double endTime, double outputTimeStep, string dir);

            double getGx();
            void setGx(double val);

            double getNU();
            void setNU(double val);

            double getH();
            void setH(double val);

            double getDY();
            void setDY(double val);

            double getR();
            void setR(double r);

            double getB();
            void setB(double val);

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
