#ifndef INVISCID_BURGERS_UPWIND_SOLVER_H
#define INVISCID_BURGERS_UPWIND_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace InviscidBurgers::Upwind {
    class Solver {
        private:
            double l;
            double dx;
            double a;
            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            void setInitialCondition(int nx, vector<double>& u);          
            
            double maxAbsDifference(int nx, vector<double>& u, vector<double>& u1);

            void updateData(int nx, vector<double>& u, vector<double>& u1);
            void writeData(vector<double>& u, double t, double dx, int nx, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

        public:
            Solver(double l, double dx, double a, double endTime, double outputTimeStep, string dir);

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
