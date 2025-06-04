#include <format>

#include <stdexcept>

#include "tridiagonal.h"

using namespace SLAE::Direct;

void calc(int n, const vector<double> &l, const vector<double> &c, const vector<double> &r, vector<double> &d, vector<double> &u) {
    if (n > 0) {
        u[0] = r[0] / c[0];
        d[0] = d[0] / c[0];

        for (int i = 1; i < n; i++) {
            auto c1 = c[i] - l[i]*u[i-1];
    
            u[i] = r[i] / c1;
            d[i] = (d[i] - l[i]*d[i-1]) / c1;
        }

        u[n-1] = d[n-1];

        for (int i = n-2; i >= 0; i--) {
            u[i] = d[i] - u[i]*u[i+1];
        }
    }
}

void Tridiagonal::solve(const vector<double> &l, const vector<double> &c, const vector<double> &r, vector<double> &d, vector<double> &u) {
    int n = check(l, c, r, d, u);

    calc(n, l, c, r, d, u);
}

int Tridiagonal::check(const vector<double> &l, const vector<double> &c, const vector<double> &r, vector<double> &d, vector<double> &u) {
    int n = c.size();

    if (l.size() != n) {
        throw runtime_error(
            format("The lengths of the l and c should be equal, but now the length of the l is {} and c is {}", l.size(), n)
        );
    }

    if (r.size() != n) {
        throw runtime_error(
            format("The lengths of the r and c should be equal, but now the length of the r is {} and c is {}", r.size(), n)
        );
    }

    if (d.size() != n) {
        throw runtime_error(
            format("The lengths of the d and c should be equal, but now the length of the d is {} and c is {}", d.size(), n)
        );
    }

    if (u.size() != n) {
        throw runtime_error(
            format("The lengths of the u and c should be equal, but now the length of the u is {} and c is {}", u.size(), n)
        );
    }

    return n;
}
