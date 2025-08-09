#ifndef POISEUILLE_FLOW_STEADY_STATE_SOR_SOLVER_H
#define POISEUILLE_FLOW_STEADY_STATE_SOR_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace PoiseuilleFlow::SteadyState::SOR {
    class Solver {
        private:
            double dx;
            double dy;

            double epsilon;

            string dir;

            path createDirectory();

            void setInitialCondition(
                int nx, int ny, 
                vector<vector<double>>& w
            );
            void setBoundaryCondition(
                int nx, int ny, 
                vector<vector<double>>& w
            );            
            
            double calculateMaxAbsResidualElement(
                int nx, int ny, 
                vector<vector<double>>& w
            );

            void writeData(
                vector<vector<double>>& w, 
                int nx, int ny, path outDir
            );

            void writeStatistics(
                vector<tuple<int, double>>& statistics, 
                path outDir
            );

        public:
            Solver(double dx, double dy, double epsilon, string dir);

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
