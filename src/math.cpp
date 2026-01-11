#include "math.h"

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

std::array<vec3d, 316> Math::getINodeDistVecs() {
    std::array<vec3d, 316> dvecs;

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

    std::sort(dirs.begin(), dirs.end(), vecLessThan);

    dirs.erase(
        std::unique(dirs.begin(), dirs.end(), vecEquals),
        dirs.end());

    return dirs;
}


void Math::buildPermutations(
    vec3d& xs, std::vector<vec3d>& permutes, 
    int idx0, int idx1) {

    if (idx0 == idx1) {
        permutes.push_back(xs);
        return;
    }

    for (int idx = idx0; idx < idx1; ++idx) {
        std::swap(xs[idx0], xs[idx]);
        buildPermutations(xs, permutes, idx0+1, idx1);
        std::swap(xs[idx0], xs[idx]);
    }

    // Remove dup permutations (if xs has dups)
    std::sort(permutes.begin(), permutes.end(), vecLessThan);

    permutes.erase(
        std::unique(permutes.begin(), permutes.end(), vecEquals),
        permutes.end());
}
