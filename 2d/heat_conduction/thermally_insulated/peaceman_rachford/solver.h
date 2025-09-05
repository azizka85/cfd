#ifndef HEAT_CONDUCTION_THERMALLY_INSULATED_PEACEMAN_RACHFORD_SOLVER_H
#define HEAT_CONDUCTION_THERMALLY_INSULATED_PEACEMAN_RACHFORD_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace HeatConduction::ThermallyInsulated::PeacemanRachford {
    class Solver {
        private:
            double alpha;

            double l;
            double h;

            double dx;
            double dy;

            double r;

            double endTime;
            double outputTimeStep;

            string dir;

            path createDirectory();

            void setInitialCondition(
                int nx, int ny, 
                vector<vector<double>>& T
            );                             
            
            double maxAbsDifference(
                int nx, int ny, 
                vector<vector<double>>& T, 
                vector<vector<double>>& T1
            );

            void updateData(
                int nx, int ny, 
                vector<vector<double>>& T, 
                vector<vector<double>>& T1
            );
            void writeData(
                vector<vector<double>>& T, 
                double t, int nx, int ny, 
                int m, path outDir
            );
            void writeStatistics(
                vector<tuple<int, double, double>>& statistics, 
                path outDir
            );

        public:
            Solver(
                double alpha, double l, double h,
                double dx, double dy, double r, 
                double endTime, double outputTimeStep, string dir
            );

            double getALPHA();
            void setALPHA(double val);

            double getL();
            void setL(double val);

            double getH();
            void setH(double val);

            double getDX();
            void setDX(double val);

            double getDY();
            void setDY(double val);

            double getR();
            void setR(double r);

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
