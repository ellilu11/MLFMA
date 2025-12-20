#pragma once

#include "math.h"

namespace Interp {
    std::pair<realVec, realVec> gaussLegendre(
        const int, const double EPS, const double, const double);

    int getNearGLNodeIdx(
        const double, const int, const double, const double);

    double evalLagrangeBasis(const double, const realVec&, const int);

    double evalTrigBasis(const double, const realVec&, const int);
}

/* [nodes, weights] = gaussLegendreTheta(l, EPS, a, b)
    * Return lth order Gauss-Legendre nodes and weights on the interval [a,b]
    * l   : quadrature order
    * EPS : minimum error to terminate Newton-Raphson
    * a   : lower bound of interval (default -1.0)
    * b   : upper bound of interval (default 1.0)
    * nodes   : Gauss-Legendre nodes
    * weights : Gauss-Legendre weights
    */
std::pair<realVec, realVec> Interp::gaussLegendre(
    const int l, const double EPS, const double a = -1.0, const double b = 1.0) {

    const double leng = b - a;
    const double mid = (a + b)/2.0;
    assert(leng > 0);

    realVec nodes(l);
    realVec weights(l);
    const int kmax = l/2; // # positive nodes = integer part of l/2

    if (l%2) { // if order is odd, middle node is at (a+b)/2
        nodes[kmax] = mid;
        auto [p, dp] = Math::legendreP(0.0, l);
        weights[kmax] = leng / (dp*dp);
    }

    for (int k = 0; k < kmax; ++k) {
        double x_k = cos(PI * (4.0*(kmax-k)-1) / (4.0*l + 2.0));
        double dp_k;
        while (true) {
            auto [p, dp] = Math::legendreP(x_k, l);
            x_k -= p/dp; // apply Newton-Raphson
            if (abs(p/dp) <= EPS) {
                dp_k = dp;
                break;
            }
        }

        const size_t kplus = l%2 ? kmax+1+k : kmax+k;
        const size_t kminus = kmax-1-k;

        nodes[kplus] = leng/2.0*x_k + mid;
        nodes[kminus] = -leng/2.0*x_k + mid;

        weights[kplus] = leng / ((1.0-x_k*x_k) * dp_k*dp_k); 
        weights[kminus] = weights[kplus];
    }

    return std::make_pair(nodes, weights);
}

/* idx = getNearGLNodeIdx(x, m)
* Get the index of the Gauss-Legendre node of order m
* nearest and less than the point x on the interval [a, b]
* x : evaluation point
* m : order of near nodes (less than order of x)
* a : lower bound of interval
* b : upper bound of interval
* idx : Index of nearest/less Gauss-Legendre node
*/
int Interp::getNearGLNodeIdx(
    const double xi, const int m, const double a = -1.0, const double b = 1.0) {

    const double leng = b - a;
    const double mid = (a + b)/2.0;
    assert(leng > 0);

    const double x = 2.0*(xi - mid) / leng; // Change to interval [-1,1]

    const int idx = m - floor(((4.0*m+2.0) * acos(x) / PI + 1.0)/4.0) - 1; // TODO: Check for rounding errors

    assert(idx >= -1 && idx < m);

    return idx;
}

/* evalLagrangeBasis(x, xs, k)
* Evaluate the Lagrange basis polynomial taking on 1 at xs[k]
* and 0 at xs[j] for j \neq k, at the point x
* x  : evaluation point
* xs : interpolation nodes
* k  : index of basis function \in {0,1,...,order}
*/
double Interp::evalLagrangeBasis(
    const double x, const realVec& xs, const int k) {

    // assert(k < xs.size());

    double product = 1.0;

    for (int j = 0; j < xs.size(); ++j) {
        if (j == k) continue;

        product *= (x - xs[j]) / (xs[k] - xs[j]);
    }

    return product;
}

/* evalTrigBasis(x, xs, k)
* Evaluate the trigonometric basis function taking on 1 at xs[k]
* and 0 at xs[j] for j \neq k, at the point x
* x  : evaluation point
* xs : interpolation nodes
* k  : index of basis function \in {0,1,...,N-1}
*/
double Interp::evalTrigBasis(
    const double x, const realVec& xs, const int k) {

    const int N = xs.size();
    assert(k <= N-1);

    const double diff = x - xs[k];
    const double denom =
        N%2 ? N*sin(diff/2.0) : N*tan(diff/2.0);

    return denom ?
        sin(N*diff/2.0) / denom :
        1.0;
}