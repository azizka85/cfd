#ifndef HEAT_CONDUCTION_THERMALLY_INSULATED_CN_TV_SOLVER_H
#define HEAT_CONDUCTION_THERMALLY_INSULATED_CN_TV_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace HeatConduction::ThermallyInsulated::CranckNicolsonWithTimeVariation {
    class Solver {
        private:
            double a0;
            double a1;
            double alpha;
            double L;
            double dx;
            double r;
            double b;
            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            void setInitialCondition(int nx, vector<double>& T);   
            
            void calculateResidualElements(int nx, double r, vector<double>& T, vector<double>& d);

            double adjustTimeStep(double t, double dt);
            
            double maxAbsDifference(int nx, vector<double>& T, vector<double>& T1);

            void updateData(int nx, vector<double>& T, vector<double>& T1);
            void writeData(vector<double>& T, double t, double dx, int nx, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

        public:
            Solver(double a0, double a1, double alpha, double L, double dx, double r, double b, double endTime, double outputTimeStep, string dir);

            double getA0();
            void setA0(double val);

            double getA1();
            void setA1(double val);

            double getAlpha();
            void setAlpha(double val);

            double getL();
            void setL(double val);

            double getDX();
            void setDX(double val);

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
