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

    // TODO: Generate recursively
    matXcd wignerD_l(const pair2d angles, const int l) {
        using namespace std;
        const auto [th, ph] = angles;

        auto sumCoeff = [th, l](int m, int n, int s) {
            int a0 = l+m-s, a1 = n-m+s, a2 = s, a3 = l-n-s;
            return pow(-1.0, n-m+s) *
                pow(cos(th/2.0), a0+a3) * pow(sin(th/2.0), a1+a2) /
                (factorial(a0) * factorial(a1) * factorial(a2) * factorial(a3));
            };

        matXcd mat = matXcd::Zero(2*l+1, 2*l+1);
        for (int n = -l; n <= l; ++n) {
            int n_ = n+l;
            double pm_n = (n < 0) ? pm(n) : 1.0;
            cmplx exp_n = expI(static_cast<double>(-n)*ph);
            for (int m = -l; m <= l; ++m) {
                int m_ = m+l;
                double pm_m = (m < 0) ? pm(m) : 1.0;
                for (int s = max(m-n, 0); s <= min(l+m, l-n); ++s)
                    mat(n_, m_) += sumCoeff(m, n, s);

                mat(n_, m_) *= exp_n
                    * pm_n / pm_m
                    * sqrt(factorial(l+n)*factorial(l-n)*factorial(l+m)*factorial(l-m));
            }
        }

        return mat;
    }
}