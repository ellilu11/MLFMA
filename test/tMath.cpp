#include <iostream>
#include "../src/interp.h"
#include "../src/interp.cpp"
#include "../src/math.cpp"
// #include "../src/fmm/node.h"

using namespace Math;

void testLegendreP(const realVec& xs, const int nmax) {
    std::cout << "\nTesting Legendre...\n";

    for (int i = 0; i < xs.size(); ++i) {
        for (int n = 0; n <= nmax; ++n) 
            std::cout << Math::legendreP(xs[i], n).first << ' ';
        
        std::cout << '\n';
    }
}

void testSphericalHankel1(const vecXd& xs, const int nmin, const int nmax) {
    std::cout << "\nTesting Spherical Hankel...\n";

    for (int i = 0; i < xs.size(); ++i) {
        for (int n = nmin; n <= nmax; ++n)
            std::cout << Math::sphericalHankel1(xs[i], n) << ' ';

        std::cout << '\n';
    }
}

std::vector<realVec> testGaussLegendre(const int maxOrder, double a, double b) {
    std::cout << "\nTesting Gauss-Legendre...\n";
    std::vector<realVec> nodesVec;

    for (int l = 1; l <= maxOrder; ++l) {
        auto [nodes, weights] = Interp::gaussLegendre(l, 1.0E-9, a, b);
        std::cout << "l = " << l << ": ";
        for (int k = 0; k < l; ++k)
            std::cout << '(' << nodes[k] << ',' << weights[k] << ") ";
        std::cout << '\n';

        nodesVec.push_back(nodes);
    }

    return nodesVec;
}

void testLagrangeInterp(const realVec& xs) {
    std::cout << "\nTesting Lagrange interpolants...\n";

    const int order = xs.size()-1;

    for (int k = 0; k <= order; ++k) {
        std::cout << "k = " << k << ": ";
        for (int j = 0; j <= order; ++j)

            std::cout << Interp::evalLagrangeBasis(xs[j], xs, k) << " ";
        std::cout << '\n';
    }
}

void testTrigInterp(const int N) {
    std::cout << "\nTesting trigonometric interpolants (N = " << N << ") ...\n";

    realVec xs;
    for (int j = 0; j < N; ++j)
        xs.push_back(2.0*PI*j/static_cast<double>(N));

    for (int k = 0; k < N; ++k) {
        std::cout << "k = " << k << ": ";
        for (int j = 0; j < N; ++j) {
            double x = xs[j];
            // if (j == k) x *= (1.0 + 1.0E-6);

            std::cout << Interp::evalTrigBasis(x, xs, k) << " ";
        }
        std::cout << '\n';
    }
}

void testNearGLNodeIdx(
    const std::vector<realVec>& nodesVec, const int m, const int lIdx, double a, double b) {

    std::cout << "\nTesting near node finding...\n";

    double x = nodesVec[m][lIdx];

    auto mIdx = Interp::getNearGLNodeIdx(x, m, a, b); // expression to test

    double y = (mIdx >= 0 ? nodesVec[m-1][mIdx] : a);

    std::cout << x << " is nearest and greater than "
        << y << ", the l = " << m << " idx = " << mIdx << " node\n";
}

void testNearPhiIdx(const double phi, const int N) {

    std::cout << "\nTesting near phi finding...\n";

    auto phi_k = [N](const int k) {
        return 2.0 * PI * k / N;
    };

    for (int i = 0; i < N; ++i)
        std::cout << phi_k(i) << ' ';
    std::cout << '\n';

    auto idx = std::floor(N * phi / (2.0 * PI)); // expression to test

    std::cout << phi << " is nearest and greater than "
        << phi_k(idx) << ", node " << idx << "\n";
}

void testINodeFuncs() {

    std::cout << "\nTesting INodeDistances()...\n";

    auto dists = Math::getINodeDistances();

    for (int i = 0; i < dists.size(); ++i)
        std::cout << i << ' ' << dists[i] << '\n';

    std::cout << "\nTesting INodeDirections()...\n";

    auto dirs = Math::getINodeDirections();

    for (int i = 0; i < dirs.size(); ++i)
        std::cout << i << ' ' << dirs[i] << '\n';

}

void testPermutation(vec3d& xs) {
    const int leng = xs.size();

    std::vector<vec3d> permutes;
    permutes.reserve(Math::factorial(leng));

    Math::buildPermutations(xs, permutes, 0, leng);
    // assert(permutes.size() == Math::factorial(leng));

    for (const auto& permute : permutes) {
        for (const auto& x : permute)
            std::cout << x << ' ';
        std::cout << '\n';
    }
}

/*void testTransl() {
    // ==================== Test translation ===================== //
    std::cout << " Testing M2L translations...\n";

    realVec testDists =
        { 2, 2.23607, 2.44949, 2.82843, 3, 3.16228, 3.31662, 3.4641,
          3.60555, 3.74166, 4.12311, 4.24264, 4.3589, 4.69042, 5.19615 };

    const auto& translTable = Node::getTables().transl;

    for (const auto& dist : testDists)
        std::cout << dist << ' ' << translTable[1].at(dist)[0] << '\n';
}*/

int main() {

    /*
    realVec xs;
    const int nmax = 10;
    for (double x = -1.0; x <= 1.0; x += 0.1 )
        xs.push_back(x);

    testLegendreP(xs, nmax);
    */

    /*
    const double nodeLeng = 5.0;
    const vecXd xs{ { 2, 2.23607, 2.44949, 2.82843, 3, 3.16228, 3.31662, 3.4641,
          3.60555, 3.74166, 4.12311, 4.24264, 4.3589, 4.69042, 5.19615 } };
    testSphericalHankel1(nodeLeng*xs, 10, 20);
    */

    /*
    const int maxOrder = 10;
    const double a = -1.0;
    const double b = 1.0;
    auto nodesVec = testGaussLegendre(maxOrder, a, b);
    */
    
    // testLagrangeInterp(nodesVec[order-1]);
    
    /*const int m = 7;
    for (int i = 0; i <=m; ++i)
        testNearGLNodeIdx(nodesVec, m, i, a, b);

    testNearPhiIdx(0, 10);*/

    vec3d xs{ 1, 1, 2 };
    testPermutation(xs);

    return 0;
}
