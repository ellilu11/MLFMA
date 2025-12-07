#pragma once

#include <cmath>
#include "types.h"

#define _USE_MATH_DEFINES

constexpr cmplx iu(0, 1);
const double PI = std::acos(-1.0);
const vec3d zeroVec = vec3d::Zero();

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

    inline size_t coord2Idx(size_t row, size_t col, int leng) {
        return row*leng + col;
    }

    inline pair2i idx2Coord(size_t idx, int leng) {
        return pair2i(idx/leng, idx%leng);
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

    inline double factorial(double n) {
        return n == 0 ? 1 : n * factorial(n-1);
    }

    inline double fallingFactorial(double x, int k) {
        return k == 0 ? 1 : x * fallingFactorial(x - 1, k - 1);
    }

    inline vec3d toSph(const vec3d& X) {
        auto x = X[0], y = X[1], z = X[2], r = X.norm();
        assert(r != 0);

        auto toPhi = [](double x, double y) {
            if (x == 0 && y == 0) return 0.0; // pick phi = 0.0
            return atan2(y, x);
            };

        return vec3d(r, std::acos(z/r), toPhi(x, y));
    }

    inline vec3d fromSph(const vec3d& R) {
        auto r = R[0], th = R[1], ph = R[2];
        return vec3d(
            r * sin(th) * cos(ph),
            r * sin(th) * sin(ph),
            r * cos(th));
    }

    inline vec3d fromCyl(const vec3d& S) {
        auto r = S[0], ph = S[1], z = S[2];
        return vec3d(
            r * cos(ph),
            r * sin(ph),
            z);
    }

    inline vec3d vecSph(const double r, const double th, const double ph) {
        return r * vec3d(
            sin(th) * cos(ph),
            sin(th) * sin(ph),
            cos(th));
    }

    inline mat3d IminusRR(const double th, const double ph) {
        auto rhat = fromSph(vec3d(1.0, th, ph));
        return mat3d::Identity() - rhat * rhat.transpose();
    }

    inline double coeffYlm(int l, int abs_m) {
        assert(abs_m <= l);
        return
            std::sqrt(factorial(l-abs_m) / static_cast<double>(factorial(l+abs_m))) * // Ylm coeffs
            pm(abs_m) * std::pow(2.0, l); // legendreLM coeffs
    }

    inline mat23d matToThPh(const double th, const double ph) {
        return mat23d{
            {  cos(th)*cos(ph),  cos(th)*sin(ph), -sin(th) },
            { -sin(ph),          cos(ph),          0.0     }
        };
    }

    inline mat3d matToSph(const double th, const double ph) {
        return mat3d{
            {  sin(th)*cos(ph),  sin(th)*sin(ph),  cos(th) },
            {  cos(th)*cos(ph),  cos(th)*sin(ph), -sin(th) },
            { -sin(ph),          cos(ph),          0.0     }
        };
    }

    inline mat3d matFromSph(const double th, const double ph) {
        return mat3d{
            {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph) },
            {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph) },
            {  cos(th),         -sin(th),          0.0     }
        };
    }

    inline mat3d rotationR(const double th, const double ph) {
        return mat3d{
            {  cos(th)*cos(ph),  cos(th)*sin(ph), -sin(th) },
            { -sin(ph),          cos(ph),          0       },
            {  sin(th)*cos(ph),  sin(th)*sin(ph),  cos(th) }
        };
    }

    inline size_t flipIdxToRange(const int i, const int size) {
        int uint_i = i;

        if (i < 0)
            uint_i = -i-1;
        else if (i >= size)
            uint_i = 2*size-i-1;

        assert(uint_i >= 0 && uint_i < size);
        return uint_i;
    }

    inline size_t wrapIdxToRange(const int i, const int size) {
        int uint_i = i;

        if (i < 0)
            uint_i = i + size;
        else if (i >= size)
            uint_i = i - size;

        assert(uint_i >= 0 && uint_i < size);
        return uint_i;
    }

    /*inline int getIdxOfVal(const double val, const realVec& list) {

    auto idx = std::find_if( list.begin(), list.end(),
        [val](double val0) {
        return std::abs(val - val0) < 1.0E-6;
        }
    );

    return idx;
    }*/

    pair2d legendreL(const double, const int);

    cmplx sphericalHankel1(const double, const int);

} // close Math::

/* legendreL(x,l)
 * Recursively evaluate the lth order Legendre polynomial
 * and its 1st derivative at the point x
 * x : evaluation point
 * l : order of Legendre polynomial
 */
pair2d Math::legendreL(const double x, const int l) {
    double p;
    double pm2 = 1.0;
    double pmm = x;

    for (int i = 2; i <= l; ++i) {
        p = ((2.0*i-1)*x*pmm - (i-1)*pm2) / static_cast<double>(i);
        pm2 = pmm;
        pmm = p;
    }

    double dp = l*(x*pmm - pm2) / (x*x - 1.0);

    return std::make_pair(p, dp);
}

/* sphericalHankel1(x,n)
     * Recursively evaluate the spherical Hankel function
     * of the 1st kind of order n at the point x
     * x : evaluation point
     * n : order of Hankel function
     */
cmplx Math::sphericalHankel1(const double x, const int n) {

    cmplx H1_nmm = -iu*expI(x) / x;

    if (!n) return H1_nmm;

    cmplx H1_n = -expI(x) * (1.0/x + iu/(x*x));

    for (int i = 0; i < n-1; ++i) {
        cmplx H1_npp = (2*n-1)/x * H1_n - H1_nmm;
        H1_nmm = H1_n;
        H1_n = H1_npp;
    }

    return H1_n;

    /*
    else if (n == 1)
        return
    else {
        assert(n >= 2);
        return
            (2*n-1)/x * sphericalHankel1(x, n-1) -
            sphericalHankel1(x, n-2);
    }*/
}