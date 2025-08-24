#ifndef POISEUILLE_FLOW_TRIANGULAR_DOMAIN_FREE_SURFACE_PEACEMAN_RACHFORD_SOLVER_H
#define POISEUILLE_FLOW_TRIANGULAR_DOMAIN_FREE_SURFACE_PEACEMAN_RACHFORD_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace PoiseuilleFlow::TriangularDomain::FreeSurface::PeacemanRachford {
    class Solver {
        private:
            double nu;
            
            double a;

            double dx;
            double dy;
            
            double r;
            
            double endTime;
            double outputTimeStep;
            
            string dir;

            tuple<double, double, double, double> generateDomain();
            tuple<bool, bool> pointInDomain(double x, double y);
            bool isPointInDomain(double x, double y);

            path createDirectory();

            void setInitialCondition(
                int nx, int ny, 
                double x0, double y0,
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
                double x0, double y0,
                int m, path outDir
            );
            void writeStatistics(
                vector<tuple<int, double, double>>& statistics, 
                path outDir
            );

        public:
            Solver(
                double nu, double a, double dx, double dy, double r, 
                double endTime, double outputTimeStep, string dir
            );

            double getNU();
            void setNU(double val);

            double getA();
            void setA(double val);

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
