#pragma once

#include <cmath>
#include <utility>
#include "types.h"

#define _USE_MATH_DEFINES

constexpr cmplx iu(0, 1);
const double PI = std::acos(-1.0);

const vec3d zeroVec = vec3d::Zero();
const std::array<vec3d, 2> poles{ vec3d(0, 0, 1), vec3d(0, 0, -1) };

namespace Math {
    constexpr double FEPS = 1.0E-6; // floating point error tolerance

    inline size_t bools2Idx(const std::array<bool,3>& x) noexcept {
        return x[0] + 2 * x[1] + 4 * x[2];
    }

    inline double sign(int k) noexcept {
        return k % 2 ? -1.0 : 1.0;
    }

    inline cmplx powI(int m) noexcept {
        switch (m % 4) {
            case 0: return 1;
            case 1: return iu;
            case 2: return -1.0;
            case 3: return -iu;
        }
    }

    inline double factorial(double n) noexcept {
        return n == 0 ? 1 : n * factorial(n-1);
    }

    inline bool approxZero(double x) noexcept {
        return fabs(x) < FEPS;
    }

    inline bool approxLess(double x, double y) noexcept {
        if (fabs(x-y) < FEPS) return false;
        return x < y;
    }

    inline bool vecEquals(const vec3d X, const vec3d Y) {
        return ((X-Y).norm()) < FEPS;
    };

    inline bool vecLessThan(const vec3d& X, const vec3d& Y) {
        if (approxLess(X[0], Y[0])) return true;
        if (approxLess(Y[0], X[0])) return false;

        if (approxLess(X[1], Y[1])) return true;
        if (approxLess(Y[1], X[1])) return false;

        return X[2] < Y[2];
    };

    inline vec3d toSph(const vec3d& X) noexcept {
        auto x = X[0], y = X[1], z = X[2], r = X.norm();
        assert(r != 0);

        auto toPhi = [](double x, double y) {
            if (x == 0 && y == 0) return 0.0; // pick phi = 0.0
            return atan2(y, x);
            };

        return vec3d(r, std::acos(z/r), toPhi(x, y));
    }

    inline vec3d fromSph(const vec3d& R) noexcept {
        auto r = R[0], th = R[1], ph = R[2];
        return r * vec3d(
            sin(th) * cos(ph),
            sin(th) * sin(ph),
            cos(th));
    }

    inline vec3d fromCyl(const vec3d& S) noexcept {
        auto r = S[0], ph = S[1], z = S[2];
        return vec3d(
            r * cos(ph),
            r * sin(ph),
            z);
    }

    inline mat23d toThPh(double th, double ph) noexcept {
        return mat23d{
            {  cos(th)*cos(ph),  cos(th)*sin(ph), -sin(th) },
            { -sin(ph),          cos(ph),          0.0     }
        };
    }

    inline mat3d ImRR(const vec3d& rhat) {
        return mat3d::Identity() - rhat * rhat.transpose();
    }

    //inline mat23d crossToThPh(double th, double ph) noexcept {
    //    return mat23d{
    //        { -sin(ph),          cos(ph),          0.0     },
    //        { -cos(th)*cos(ph), -cos(th)*sin(ph),  sin(th) }
    //    };
    //}

    //inline mat32d matFromThPh(double th, double ph) noexcept {
    //    return mat32d{
    //        {  cos(th)*cos(ph), -sin(ph) },
    //        {  cos(th)*sin(ph),  cos(ph) },
    //        { -sin(th),          0.0     }
    //    };
    //}

    //inline mat3d matToSph(double th, double ph) noexcept {
    //    return mat3d{
    //        {  sin(th)*cos(ph),  sin(th)*sin(ph),  cos(th) },
    //        {  cos(th)*cos(ph),  cos(th)*sin(ph), -sin(th) },
    //        { -sin(ph),          cos(ph),          0.0     }
    //    };
    //}

    //inline mat3d matFromSph(double th, double ph) noexcept {
    //    return mat3d{
    //        {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph) },
    //        {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph) },
    //        {  cos(th),         -sin(th),          0.0     }
    //    };
    //}

    inline Eigen::Matrix3cd dyadicG(const vec3d& X, double k) noexcept {
        const double r = X.norm(), kr = k*r, invkrsq = 1.0/(kr*kr);
        const cmplx iinvkr = iu/kr;
        const vec3d& rhat = X / r; // X.normalized()
        const mat3d& RR = rhat * rhat.transpose();

        return
            exp(iu*kr)/r * (
                mat3d::Identity() * (1.0 + iinvkr - invkrsq) -
                RR * (1.0 + 3.0*iinvkr - 3.0*invkrsq));
    };

    inline size_t flipIdxToRange(int i, int size) noexcept {
        int uint_i = i;

        if (i < 0)
            uint_i = -i-1;
        else if (i >= size)
            uint_i = 2*size-i-1;

        assert(uint_i >= 0 && uint_i < size);
        return uint_i;
    }

    inline size_t wrapIdxToRange(int i, int size) noexcept {
        int uint_i = i;

        if (i < 0)
            uint_i = i + size;
        else if (i >= size)
            uint_i = i - size;

        assert(uint_i >= 0 && uint_i < size);
        return uint_i;
    }

    inline pair2cd givensRotation(cmplx z, cmplx w) noexcept {
        // const double norm = sqrt(std::norm(z) + std::norm(w)); // squared norms
        const cmplx t = sqrt(z*z + w*w);
        return make_pair(z/t, w/t);
    }

    vec3d idx2pm(int);

    pair2d legendreP(double, int);

    cmplx sphericalHankel1(double, int);

    realVec getINodeDistances();

    std::array<vec3d, 316> getINodeDistVecs();

    std::vector<vec3d> getINodeDirections();

    void buildPermutations(vec3d&, std::vector<vec3d>&, int, int);
}


