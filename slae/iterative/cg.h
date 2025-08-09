#ifndef SLAE_ITERATIVE_CG_H
#define SLAE_ITERATIVE_CG_H

#include <vector>

#include <functional>

using namespace std;

namespace SLAE::Iterative::CG {
    void preconditionedCG(
        function<
            void (
                const vector<double>& r, 
                vector<double>& y
            )
        > preconditioner,
        function<void (int iter, double maxAbsRes)> statistics,
        const vector<int>& ip, const vector<int>& jp,
        const vector<double>& e, vector<double> &b,
        vector<double>& x, double epsilon
    );
}

#endif
