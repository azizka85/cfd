#ifndef SLAE_UTILS_H
#define SLAE_UTILS_H

#define VAR_NAME(x) #x

#include <vector>

using namespace std;

namespace SLAE::Utils {
    void check(int n, int m, const char* name);

    double dot(const vector<double>& a, const vector<double>& b);

    void csrMultVec(
        const vector<int>& ip, const vector<int>& jp, 
        const vector<double>& e, const vector<double>& v,
        vector<double>& r
    );

    void csrSymMultVec(
        const vector<int>& ip, const vector<int>& jp, 
        const vector<double>& e, const vector<double>& v,
        vector<double>& r
    );

    void axpby(
        double a, double b, 
        const vector<double>& x, const vector<double>& y,
        vector<double>& r
    );

    void update(double a, const vector<double>& x, vector<double>& y);

    double maxAbsElem(const vector<double>& a);
}

#endif
