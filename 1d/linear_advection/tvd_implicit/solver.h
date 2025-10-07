#ifndef LINEAR_ADVECTION_TVD_IMPLICIT_SOLVER_H
#define LINEAR_ADVECTION_TVD_IMPLICIT_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace LinearAdvection::TVD::Implicit {
    class Solver {
        private:
            double u;
            double d;
            double l;
            double dx;
            double a;
            double epsilon;
            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            void setInitialCondition(int nx, vector<double>& u);          
            
            double maxAbsDifference(int nx, vector<double>& u, vector<double>& u1);

            void updateData(int nx, vector<double>& u, vector<double>& u1);
            void writeData(vector<double>& u, double t, double dx, int nx, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

            double q(double x);

        public:
            Solver(double u, double d, double l, double dx, double a, double epsilon, double endTime, double outputTimeStep, string dir);

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
