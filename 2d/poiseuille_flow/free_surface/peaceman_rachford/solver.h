#ifndef POISEUILLE_FLOW_FREE_SURFACE_PEACEMAN_RACHFORD_SOLVER_H
#define POISEUILLE_FLOW_FREE_SURFACE_PEACEMAN_RACHFORD_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace PoiseuilleFlow::FreeSurface::PeacemanRachford {
    class Solver {
        private:
            double nu;
            double dx;
            double dy;
            double r;
            double endTime;
            double outputTimeStep;
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
            
            double maxAbsDifference(
                int nx, int ny, 
                vector<vector<double>>& w, 
                vector<vector<double>>& w1
            );

            void updateData(
                int nx, int ny, 
                vector<vector<double>>& w, 
                vector<vector<double>>& w1
            );
            void writeData(
                vector<vector<double>>& w, 
                double t, int nx, int ny, 
                int m, path outDir
            );
            void writeStatistics(
                vector<tuple<int, double, double>>& statistics, 
                path outDir
            );

        public:
            Solver(
                double nu, double dx, double dy, double r, 
                double endTime, double outputTimeStep, string dir
            );

            double getNU();
            void setNU(double val);

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
