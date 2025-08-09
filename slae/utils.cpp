#include <format>

#include <stdexcept>

#include "utils.h"

using namespace SLAE;

void Utils::check(int n, int m, const char* name) {
    if (m != n) {
        throw runtime_error(
            format(
                "The lengths of the vector {} should be {}, but now the length is {}", 
                name, n, m
            )
        );
    }
}

double Utils::dot(const vector<double> &a, const vector<double> &b) {
    int n = a.size();

    check(n, b.size(), VAR_NAME(b));

    double r = 0;

    for (int i = 0; i < n; i++) {
        r += a[i]*b[i];
    }

    return r;
}

void Utils::csrMultVec(
    const vector<int> &ip, const vector<int> &jp, 
    const vector<double> &e, const vector<double> &v,
    vector<double>& r
) {
    int n = ip.size();

    if (n > 0) {
        if (ip[0] != 0) {
            throw runtime_error(
                format(
                    "The first element of the vector {} should be 0, but now it is {}", 
                    VAR_NAME(ip), ip[0]
                )
            );
        }

        check(n - 1, v.size(), VAR_NAME(v));
        check(n - 1, r.size(), VAR_NAME(r));
        check(ip[n-1], jp.size(), VAR_NAME(jp));
        check(ip[n-1], e.size(), VAR_NAME(e));

        for (int i = 1; i < n; i++) {
            r[i-1] = 0;

            for (int j = ip[i-1]; j < ip[i]; j++) {
                r[i-1] += e[j]*v[jp[j]];    
            }
        }
    }
}

void SLAE::Utils::csrSymMultVec(
    const vector<int> &ip, const vector<int> &jp, 
    const vector<double> &e, const vector<double> &v, 
    vector<double> &r
) {
    int n = ip.size();

    if (n > 0) {
        if (ip[0] != 0) {
            throw runtime_error(
                format(
                    "The first element of the vector {} should be 0, but now it is {}", 
                    VAR_NAME(ip), ip[0]
                )
            );
        }

        check(n - 1, v.size(), VAR_NAME(v));
        check(n - 1, r.size(), VAR_NAME(r));
        check(ip[n-1], jp.size(), VAR_NAME(jp));
        check(ip[n-1], e.size(), VAR_NAME(e));

        for (int i = 1; i < n; i++) {
            r[i-1] = 0;

            for (int j = ip[i-1]; j < ip[i]; j++) {
                int k = jp[j];

                r[i-1] += e[j]*v[k];    

                if (k != i - 1) {
                    r[k] += e[j]*v[i-1];
                }
            }
        }
    }
}

void Utils::axpby(
    double a, double b, 
    const vector<double> &x, const vector<double> &y, 
    vector<double> &r
) {
    int n = x.size();

    check(n, y.size(), VAR_NAME(y));
    check(n, r.size(), VAR_NAME(r));

    for (int i = 0; i < n; i++) {
        r[i] = a*x[i] + b*y[i];
    }
}

void Utils::update(double a, const vector<double> &x, vector<double> &y) {
    int n = x.size();

    check(n, y.size(), VAR_NAME(y));

    for (int i = 0; i < n; i++) {
        y[i] = a*x[i];
    }
}

double Utils::maxAbsElem(const vector<double> &a) {
    double maxElem = 0;

    for (auto e: a) {
        maxElem = max(maxElem, abs(e));
    }

    return maxElem;
}
