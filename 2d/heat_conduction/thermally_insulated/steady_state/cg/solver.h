#ifndef HEAT_CONDUCTION_THERMALLY_INSULATED_STEADY_STATE_CG_SOLVER_H
#define HEAT_CONDUCTION_THERMALLY_INSULATED_STEADY_STATE_CG_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace HeatConduction::ThermallyInsulated::SteadyState::CG {
    class Solver {
        private:
            double l;
            double h;

            double dx;
            double dy;

            double epsilon;

            string dir;

            path createDirectory();

            void setInitialCondition(int nx, int ny, vector<double>& T);
            
            tuple<vector<int>, vector<int>, vector<double>> buildMatrix(int nx, int ny);
            void calculateResidualElements(int nx, int ny, vector<double>& T, vector<double>& d);

            void writeData(vector<double>& T, int nx, int ny, path outDir);

            void writeStatistics(
                vector<tuple<int, double>>& statistics, 
                path outDir
            );

        public:
            Solver(
                double l, double h,
                double dx, double dy, double epsilon, 
                string dir
            );

            double getL();
            void setL(double val);

            double getH();
            void setH(double val);

            double getDX();
            void setDX(double val);

            double getDY();
            void setDY(double val);

            double getEpsilon();
            void setEpsilon(double val);

            string getDir();                                                                        
            void setDir(string val);

            void solve();
    };  
}

#endif
