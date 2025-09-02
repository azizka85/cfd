#ifndef POISEUILLE_FLOW_STEADY_STATE_ELLIPTIC_DOMAIN_RICHARDSON_SOLVER_H
#define POISEUILLE_FLOW_STEADY_STATE_ELLIPTIC_DOMAIN_RICHARDSON_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace PoiseuilleFlow::SteadyState::EllipticDomain::Richardson {
    class Solver {
        private:
            double g;
            double nu;
            
            double a;
            double b;

            double dx;
            double dy;

            double epsilon;

            string dir;

            tuple<double, double, double, double> generateDomain();
            tuple<bool, bool> pointInDomain(double x, double y);

            path createDirectory();

            void setInitialCondition(
                int nx, int ny, 
                vector<vector<double>>& w
            );           
            
            double calculateMaxAbsResidualElement(
                int nx, int ny, 
                double x0, double y0,
                vector<vector<double>>& w
            );

            void writeData(
                vector<vector<double>>& w, 
                int nx, int ny, 
                double x0, double y0,
                path outDir
            );

            void writeStatistics(
                vector<tuple<int, double>>& statistics, 
                path outDir
            );

        public:
            Solver(
                double g, double nu, 
                double a, double b,
                double dx, double dy, double epsilon, 
                string dir
            );

            double getG();
            void setG(double val);

            double getNU();
            void setNU(double val);

            double getA();
            void setA(double val);

            double getB();
            void setB(double val);

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
