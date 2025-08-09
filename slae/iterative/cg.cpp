#include <slae/utils.h>

#include "cg.h"

using namespace SLAE::Utils;
using namespace SLAE::Iterative;

void CG::preconditionedCG(
    function<
        void(
            const vector<double> &r, 
            vector<double> &y
        )
    > preconditioner, 
    function<void(int iter, double maxAbsRes)> statistics, 
    const vector<int> &ip, const vector<int> &jp, 
    const vector<double> &e, vector<double> &b,
    vector<double> &x, double epsilon
) {
    int n = x.size();

    vector<double> p(n);
    vector<double> r(n);
    
    int iter = 1;

    csrSymMultVec(ip, jp, e, x, r);
    axpby(1, -1, r, b, r);

    if (preconditioner) {
        preconditioner(r, b);
    } else {
        update(1, r, b);
    }

    update(-1, b, p);

    double rb = dot(r, b);

    csrSymMultVec(ip, jp, e, p, b);

    double pb = dot(p, b);
    double error = maxAbsElem(r);

    while (error > epsilon) {
        double alpha = rb / pb;

        axpby(1, alpha, x, p, x);
        axpby(1, alpha, r, b, r);

        if (preconditioner) {
            preconditioner(r, b);
        } else {
            update(1, r, b);
        }

        double rb1 = dot(r, b);
        double beta = rb1 / rb;

        rb = rb1;

        axpby(-1, beta, b, p, p);
        csrSymMultVec(ip, jp, e, p, b);

        pb = dot(p, b);
        error = maxAbsElem(r);

        if (statistics) {
            statistics(iter, error);
        }

        iter += 1;
    }
}
