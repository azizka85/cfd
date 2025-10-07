#ifndef LINEAR_ADVECTION_ENO_X2T2_SOLVER_H
#define LINEAR_ADVECTION_ENO_X2T2_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace LinearAdvection::ENO::X2T2 {
    class Solver {
        private:
            double u;
            double d;
            double l;
            double dx;
            double a;
            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            void setInitialCondition(int nx, vector<double>& C);        
            
            double maxAbsDifference(int nx, vector<double>& C, vector<double>& C1);

            void updateData(int nx, vector<double>& C, vector<double>& C1);
            void writeData(vector<double>& C, double t, double dx, int nx, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

        public:
            Solver(double u, double d, double l, double dx, double a, double endTime, double outputTimeStep, string dir);

            double getU();
            void setU(double val);

            double getD();
            void setD(double val);

            double getL();
            void setL(double val);

            double getDX();
            void setDX(double val);

            double getA();
            void setA(double val);

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
