#pragma once

#include <cmath>
#include "types.h"

#define _USE_MATH_DEFINES

constexpr cmplx iu(0, 1);
const double PI = std::acos(-1.0);

const vec3d zeroVec = vec3d::Zero();
const vec3d northVec(0, 0, 1);
const vec3d southVec(0, 0, -1);

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

    inline bool approxZero(double x) noexcept {
        return fabs(x) < FEPS;
    }

    inline bool approxLess(double x, double y) noexcept {
        if (fabs(x-y) < FEPS) return false;
        return x < y;
    }

    //inline bool vecEquals(const vec3d X, const vec3d Y) noexcept {
    //    return ((X-Y).norm()) < 3.0*FEPS;
    //}

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

} // namespace Math

/* idx2pm(x)
 * Convert x to reverse binary, then replace each bit as 0 -> -1, 1 -> 1
 * returning the result
 * x : integer \in {0, ... , 7}
 */
vec3d Math::idx2pm(int x) {
    assert(x < 8);
    auto xmod4 = x%4;
    vec3d bits(xmod4%2, xmod4/2, x/4);

    for (auto& bit : bits)
        bit = (bit == 0 ? -1 : 1);

    return bits;
}

/* legendreL(x,l)
 * Recursively evaluate the lth order Legendre polynomial
 * and its 1st derivative at the point x
 * x : evaluation point
 * l : order of Legendre polynomial
 */
pair2d Math::legendreP(double x, int n) {

    double P_nmm = 1.0;

    if (!n) return std::make_pair(P_nmm, 0.0);

    double P_n = x;

    for (int i = 1; i < n; ++i) {
        double P_npp = ((2.0*i+1)*x*P_n - i*P_nmm) / static_cast<double>(i+1);
        P_nmm = P_n;
        P_n = P_npp;
    }

    double dP_n = n*(x*P_n - P_nmm) / (x*x - 1.0);

    return std::make_pair(P_n, dP_n);
}

/*pair2d Math::legendreL(double x, int l) {
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
}*/

/* sphericalHankel1(x,n)
 * Recursively evaluate the spherical Hankel function
 * of the 1st kind of order n at the point x
 * x : evaluation point
 * n : order of Hankel function
 */
cmplx Math::sphericalHankel1(double x, int n) { // TODO: Double check

    cmplx H1_nmm = -iu*exp(iu*x) / x;

    if (!n) return H1_nmm;

    cmplx H1_n = -exp(iu*x) * (1.0/x + iu/(x*x));

    for (int i = 1; i < n; ++i) {
        cmplx H1_npp = (2.0*i + 1.0)/x * H1_n - H1_nmm;
        H1_nmm = H1_n;
        H1_n = H1_npp;
    }

    return H1_n;
}

realVec Math::getINodeDistances() {
    realVec dists;

    for (double dz = 0; dz <= 3; ++dz)
        for (double dy = 0; dy <= 3; ++dy)
            for (double dx = 0; dx <= 3; ++dx) {
                auto dir = vec3d(dx, dy, dz);
                auto dist = dir.norm();

                if (dir.lpNorm<Eigen::Infinity>() > 1.0)
                    dists.push_back(dist);

            }

    std::sort(dists.begin(), dists.end());

    dists.erase(std::unique(dists.begin(), dists.end()), dists.end());

    // for (const auto& dist : dists) std::cout << dist << "\n";

    return dists;
}

std::array<vec3d,316> Math::getINodeDistVecs() {
    std::array<vec3d,316> dvecs;

    int idx = 0;
    for (double dz = -3; dz <= 3; ++dz)
        for (double dy = -3; dy <= 3; ++dy)
            for (double dx = -3; dx <= 3; ++dx) {
                auto dvec = vec3d(dx, dy, dz);

                if (dvec.lpNorm<Eigen::Infinity>() > 1.0)
                    dvecs[idx++] = dvec;
            }

    // for (const auto& dvec : dvecs) std::cout << dvec.norm() << "\n";

    return dvecs;
}

std::vector<vec3d> Math::getINodeDirections() {
    std::vector<vec3d> dirs(316);

    int idx = 0;
    for (double dz = -3; dz <= 3; ++dz)
        for (double dy = -3; dy <= 3; ++dy)
            for (double dx = -3; dx <= 3; ++dx) {
                auto dir = vec3d(dx, dy, dz);
                auto dist = dir.norm();

                if (dir.lpNorm<Eigen::Infinity>() > 1.0) 
                    dirs[idx++] = dir/dist;
            }

    auto vecLessThan = [&](const vec3d& X, const vec3d& Y) {
        if (approxLess(X[0], Y[0])) return true;
        if (approxLess(Y[0], X[0])) return false;

        if (approxLess(X[1], Y[1])) return true;
        if (approxLess(Y[1], X[1])) return false;

        return X[2] < Y[2];
        };

    auto vecEquals = [](const vec3d X, const vec3d Y) {
        return ((X-Y).norm()) < FEPS;
        };

    std::sort(dirs.begin(), dirs.end(), vecLessThan);

    dirs.erase(
        std::unique(dirs.begin(), dirs.end(), vecEquals),
        dirs.end());

    return dirs;
}