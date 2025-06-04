#ifndef SLAE_DIRECT_TRIDIAGONAL_H
#define SLAE_DIRECT_TRIDIAGONAL_H

#include <vector>

using namespace std;

namespace SLAE::Direct::Tridiagonal {
    void solve(const vector<double>& l, const vector<double>& c, const vector<double>& r, vector<double>& d, vector<double>& u);

    int check(const vector<double>& l, const vector<double>& c, const vector<double>& r, vector<double>& d, vector<double>& u);
}

#endif
