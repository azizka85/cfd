#include <format>

#include <stdexcept>

#include "augmented_tridiagonal.h"

using namespace SLAE::Direct;

void calc(
    int n, 
    const vector<double> &l, const vector<double> &c, const vector<double> &r, 
    const vector<double> &e, const vector<double> &f, 
    double g, double h,
    vector<double>& e1, vector<double> &d, vector<double> &T
) {
    if (n > 0) {
        double f1 = f[0];

        T[0] = r[0] / c[0];
        d[0] = d[0] / c[0];
        e1[0] = e[0] / c[0];

        if (n > 1) {
            for (int i = 1; i < n-1; i++) {
                auto c1 = c[i] - l[i]*T[i-1];

                T[i] = r[i] / c1;
                d[i] = (d[i] - l[i]*d[i-1]) / c1;
                e1[i] = (e[i] - l[i]*e1[i-1]) / c1;

                g -= f1 * e1[i-1];
                h -= f1 * d[i-1];
                f1 = f[i] - f1*T[i-1];            
            }

            auto c1 = c[n-1] - l[n-1]*T[n-2];

            g -= f1 * e1[n-2];
            h -= f1 * d[n-2];
            f1 = f[n-1] - f1*T[n-2];

            double lambda = h / g;

            if (abs(c1) < 1e-12) {
                lambda = (d[n-1] - l[n-1]*d[n-2]) / (e[n-1] - l[n-1]*e1[n-2]);

                g /= f1;
                h /= f1;

                T[n-1] = h - lambda * g;
            } else {
                d[n-1] = (d[n-1] - l[n-1]*d[n-2]) / c1;
                e1[n-1] = (e[n-1] - l[n-1]*e1[n-2]) / c1;

                g -= f1 * e1[n-1];
                h -= f1 * d[n-1];

                lambda = h / g;

                T[n-1] = d[n-1] - lambda * e1[n-1];
            }

            for (int i = n-2; i >= 0; i--) {
                T[i] = d[i] - T[i] * T[i+1] - lambda * e1[i];
            }
        }
    }
}

void AugmentedTridiagonal::solve(
    const vector<double> &l, const vector<double> &c, const vector<double> &r, 
    const vector<double> &e, const vector<double> &f, 
    double g, double h,
    vector<double>& e1, vector<double> &d, vector<double> &T
) {
    int n = check(l, c, r, e, f, e1, d, T);

    calc(n, l, c, r, e, f, g, h, e1, d, T);
}

int AugmentedTridiagonal::check(
    const vector<double> &l, const vector<double> &c, const vector<double> &r, 
    const vector<double> &e, const vector<double> &f, 
    vector<double>& e1, vector<double> &d, vector<double> &T
) {
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

    if (e.size() != n) {
        throw runtime_error(
            format("The lengths of the e and c should be equal, but now the length of the e is {} and c is {}", e.size(), n)
        );
    } 

    if (f.size() != n) {
        throw runtime_error(
            format("The lengths of the f and c should be equal, but now the length of the f is {} and c is {}", f.size(), n)
        );
    } 

    if (e1.size() != n) {
        throw runtime_error(
            format("The lengths of the g and c should be equal, but now the length of the g is {} and c is {}", e1.size(), n)
        );
    } 

    if (d.size() != n) {
        throw runtime_error(
            format("The lengths of the d and c should be equal, but now the length of the d is {} and c is {}", d.size(), n)
        );
    }    

    if (T.size() != n) {
        throw runtime_error(
            format("The lengths of the u and c should be equal, but now the length of the u is {} and c is {}", T.size(), n)
        );
    }

    return n;
}

void AugmentedTridiagonal::solve(
    const double l, 
    const double c, const double c0, 
    const double r, const double r0, 
    const double f, const double f0, const double f1, 
    const double h, vector<double> &T
) {
    int n = T.size();
    
    if (n > 0) {
        double ft = f0;

        T[0] = r0 / c0;

        if (n > 1) {
            for (int i = 1; i < n-1; i++) {
                auto ct = c - l*T[i-1];

                T[i] = r / ct;
                
                ft = f - ft*T[i-1];            
            }            

            ft = f1 - ft*T[n-2];

            T[n-1] = h / ft;
            
            for (int i = n-2; i >= 0; i--) {
                T[i] = -T[i] * T[i+1];
            }
        }
    }
}
