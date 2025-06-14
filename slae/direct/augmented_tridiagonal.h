#ifndef SLAE_DIRECT_AUGMENTED_TRIDIAGONAL_H
#define SLAE_DIRECT_AUGMENTED_TRIDIAGONAL_H

#include <vector>

using namespace std;

namespace SLAE::Direct::AugmentedTridiagonal {
    void solve(
        const vector<double>& l, const vector<double>& c, const vector<double>& r, 
        const vector<double>& e, const vector<double>& f, 
        double g, double h,
        vector<double>& e1, vector<double>& d, vector<double>& T
    );   

    int check(
        const vector<double>& l, const vector<double>& c, const vector<double>& r, 
        const vector<double>& e, const vector<double>& f,         
        vector<double>& e1, vector<double>& d, vector<double>& T
    );

    void solve(
        const double l, 
        const double c, const double c0, 
        const double r, const double r0, 
        const double f, const double f0, const double f1, 
        const double h, vector<double> &T
    );
}

#endif
