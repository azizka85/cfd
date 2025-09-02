#ifndef POISEUILLE_FLOW_STEADY_STATE_TRIANGULAR_DOMAIN_FREE_SURFACE_CG_SOLVER_H
#define POISEUILLE_FLOW_STEADY_STATE_TRIANGULAR_DOMAIN_FREE_SURFACE_CG_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace PoiseuilleFlow::SteadyState::TriangularDomain::FreeSurface::CG {
    class Solver {
        private:
            double a;

            double dx;
            double dy;

            double epsilon;

            string dir;

            tuple<double, double, double, double> generateDomain();
            tuple<bool, bool> pointInDomain(double x, double y);
            bool isPointInDomain(double x, double y);

            path createDirectory();

            void setInitialCondition(int nx, int ny,  vector<double>& w);               

            tuple<vector<int>, vector<int>, vector<double>> buildMatrix(int nx, int ny, double x0, double y0);
            void calculateResidualElements(int nx, int ny, double x0, double y0, vector<double>& w, vector<double>& d);

            void writeData(vector<double>& w, int nx, int ny, double x0, double y0, path outDir);
            void writeStatistics(vector<tuple<int, double>>& statistics, path outDir);

        public:
            Solver(
                double a, double dx, double dy, double epsilon,
                string dir
            );

            double getA();
            void setA(double val);

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
