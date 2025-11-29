#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES

using realVec = std::vector<double>;
using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;

using pair2i = std::pair<int, int>;
using pair2d = std::pair<double, double>;

using vec3i = Eigen::Vector3i;
using vec3d = Eigen::Vector3d;
using vec3cd = Eigen::Vector3cd;
using vecXcd = Eigen::VectorXcd;

using mat3d = Eigen::Matrix3d;
using matXcd = Eigen::MatrixXcd;

constexpr cmplx iu(0, 1);
const double PI = std::acos(-1.0);
const vec3d zeroVec = vec3d::Zero();

template <typename T>
std::vector<T> operator+ (const std::vector<T>& zs, const std::vector<T>& ws) {
    std::vector<T> sum;
    for (size_t i = 0; i < zs.size(); ++i)
        sum.push_back(zs[i] + ws[i]);
    return sum;
}

std::ostream& operator<< (std::ostream& os, const vec3d& X) {
    os << X[0] << " " << X[1] << " " << X[2];
    return os;
}

std::istream& operator>>(std::istream& is, vec3d& X) {
    double x, y, z;
    if (is >> x >> y >> z)
        X = vec3d{ x, y, z };
    return is;
}

std::istream& operator>>(std::istream& is, vec3i& X) {
    int x, y, z;
    if (is >> x >> y >> z)
        X = vec3i{ x, y, z };
    return is;
}

std::istream& operator>>(std::istream& is, Eigen::Vector4i& X) {
    int x, y, z, w;
    if (is >> x >> y >> z >> w)
        X = Eigen::Vector4i{ x, y, z, w };
    return is;
}

/*
template <typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& X) {
    T x, y, z;
    if (is >> x >> y >> z)
        X = std::vector<T>{ x, y, z };
    return is;
}*/

std::array<bool, 3> operator> (const vec3d& x, const vec3d& y) {
    std::array<bool, 3> bools{ x[0] > y[0], x[1] > y[1], x[2] > y[2] };
    return bools;
}

namespace Math {

    /* bools2Idx(x)
     * Convert bools into branchIdx \in {0, ... ,7}
     */
    inline size_t bools2Idx(const std::array<bool, 3>& x) {
        return x[0] + 2 * x[1] + 4 * x[2];
    }

    /* idx2pm(x)
     * Convert x to reverse binary, then replace each bit as 0 -> -1, 1 -> 1
     * returning the result
     * x : integer \in {0, ... , 7}
     */
    inline vec3d idx2pm(const int x) {
        assert(x < 8);
        auto xmod4 = x%4;
        vec3d bits(xmod4%2, xmod4/2, x/4);

        for (auto& bit : bits)
            bit = (bit == 0 ? -1 : 1);

        return bits;
    }

    inline double pm(const int k) {
        return k % 2 ? -1.0 : 1.0;
    }

    inline cmplx expI(const double arg) {
        return std::exp(iu*arg);
    }

    inline cmplx powI(const uint32_t m) {
        switch (m % 4) {
            case 0: return 1;
            case 1: return iu;
            case 2: return -1.0;
            case 3: return -iu;
        }
    }

    //const uint64_t factorial(int n) {
    //    return n == 0 ? 1 : n * factorial(n-1);
    //}

    double factorial(double n) {
        return n == 0 ? 1 : n * factorial(n-1);
    }

    double fallingFactorial(double x, int k) {
        return k == 0 ? 1 : x * fallingFactorial(x - 1, k - 1);
    }

    inline vec3d fromSph(const vec3d& R) {
        auto r = R[0], th = R[1], ph = R[2];
        return vec3d(
            r * std::sin(th) * std::cos(ph),
            r * std::sin(th) * std::sin(ph),
            r * std::cos(th));
    }

    inline vec3d fromCyl(const vec3d& S) {
        auto r = S[0], ph = S[1], z = S[2];
        return vec3d(
            r * std::cos(ph),
            r * std::sin(ph),
            z);
    }

    inline vec3d toSph(const vec3d& X) {
        auto x = X[0], y = X[1], z = X[2], r = X.norm();
        assert(r != 0);

        auto toPhi = [](double x, double y) {
            if (x == 0 && y == 0) return 0.0; // pick phi = 0.0
            return std::atan2(y, x);
            };

        return vec3d(r, std::acos(z/r), toPhi(x, y));
    }

    inline double coeffYlm(int l, int abs_m) {
        assert(abs_m <= l);
        return
            std::sqrt(factorial(l-abs_m) / static_cast<double>(factorial(l+abs_m))) * // Ylm coeffs
            pm(abs_m) * std::pow(2.0, l); // legendreLM coeffs
    }

    inline mat3d matFromSph(const double th, const double ph) {
        return mat3d{
            {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph)/sin(th) },
            {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph)/sin(th) },
            {  cos(th),         -sin(th),          0.0             }
        };
    }

    inline mat3d rotationR(const double th, const double ph) {
        return mat3d{
            {  cos(th)*cos(ph),  cos(th)*sin(ph), -sin(th) },
            { -sin(ph),          cos(ph),          0       },
            {  sin(th)*cos(ph),  sin(th)*sin(ph),  cos(th) }
        };
    }

    /* legendreL(x,l)
     * Recursively evaluate the lth order Legendre polynomial 
     * and its 1st derivative at the point x
     * x : evaluation point
     * l : order of Legendre polynomial
     */
    pair2d legendreL(double x, int l) {
        double p;
        double pm2 = 1.0;
        double pmm = x;

        for (int i = 2; i <= l; ++i) {
            p = ((2.0*i-1)*x*pmm - (i-1)*pm2) / static_cast<double>(i);
            pm2 = pmm;
            pmm = p;
        }

        double dp = l*(x*pmm - pm2) / (x*x - 1.0);

        return pair2d(p, dp);
    }

    /* gaussLegendreTheta(l)
    * Return lth order Gauss-Legendre nodes and weights on the interval [0,\pi]
    * l : quadrature order
    * EPS : minimum error to terminate Newton-Raphson
    */
    std::pair<realVec,realVec> gaussLegendreTheta(int l, const double EPS) {
        realVec nodes(l);
        realVec weights(l);
        const int kmax = l/2; // # positive nodes = integer part of l/2

        if (l%2) { // if order is odd, middle node is at \pi/2
            nodes[kmax] = PI / 2.0;
            auto [p, dp] = legendreL(0.0, l);
            weights[kmax] = PI / (dp*dp);
        }

        for (int k = 0; k < kmax; ++k) {
            double x_k = cos(PI * (4.0*(kmax-k)-1) / (4.0*l + 2.0));
            double dp_k;
            while (1) {
                auto [p, dp] = legendreL(x_k, l);
                x_k -= p/dp; // apply Newton-Raphson
                if (abs(p/dp) <= EPS) {
                    dp_k = dp;
                    break;
                }
            }

            const size_t kplus = l%2 ? kmax+1+k : kmax+k;
            const size_t kminus = kmax-1-k;
                
            nodes[kplus] = PI/2.0*(x_k + 1.0);
            nodes[kminus] = PI/2.0*(-x_k + 1.0);

            weights[kplus] = PI / ((1.0-x_k*x_k) * dp_k*dp_k);
            weights[kminus] = weights[kplus];
        }

        return std::pair<realVec,realVec>(nodes, weights);
    }

} // close Math::