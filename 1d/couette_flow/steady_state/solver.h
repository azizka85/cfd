#ifndef COUETTE_FLOW_STEADY_STATE_SOLVER_H
#define COUETTE_FLOW_STEADY_STATE_SOLVER_H

#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace CouetteFlow::SteadyState {
    class Solver {
        private:
            double uTop;
            double h;
            double dy;
            string dir;

            path createDirectory();

            void setBoundaryCondition(int ny, vector<double>& u);

            void writeData(vector<double>& u, double dy, int ny, path outDir);

        public:
            Solver(double uTop, double h, double dy, string dir);

            double getUTop();
            void setUTop(double val);

            double getH();
            void setH(double val);

            double getDY();
            void setDY(double val);

            string getDir();                                                                        
            void setDir(string val);

            void solve();
    };
}

#endif
