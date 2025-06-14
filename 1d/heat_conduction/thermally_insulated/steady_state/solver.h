#ifndef HEAT_CONDUCTION_THERMALLY_INSULATED_1D_STEADY_STATE_SOLVER_H
#define HEAT_CONDUCTION_THERMALLY_INSULATED_1D_STEADY_STATE_SOLVER_H

#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace HeatConduction::ThermallyInsulated::SteadyState {
    class Solver {
        private:
            double a;
            double L;
            double dx;
            string dir;

            path createDirectory();

            void writeData(vector<double>& T, double dx, int nx, path outDir);

        public:
            Solver(double a, double L, double dx, string dir);

            double getA();
            void setA(double val);

            double getL();
            void setL(double val);

            double getDX();
            void setDX(double val);

            string getDir();                                                                        
            void setDir(string val);

            void solve();
    };
}

#endif
