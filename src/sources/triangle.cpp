#include "triangle.h"

using namespace std;

Triangle::Triangle(
    const vec3i& vIdx,
    const std::vector<vec3d>& vList,
    const Precision quadPrec)
    : vIdx(vIdx),
    Xs({ vList[vIdx[0]], vList[vIdx[1]], vList[vIdx[2]] }),
    center( (Xs[0]+Xs[1]+Xs[2])/3.0 )
{
    for (int i = 0; i < 3; ++i) {
        const int ipp = Math::wrapIdxToRange(i+1, 3);
        Ds[i] = Xs[ipp] - Xs[i];
    }

    nhat = (Ds[0].cross(Ds[1])).normalized();

    buildQuads(quadPrec);
};

int Triangle::prec2Int(Precision quadPrec) {
    return
        [&]() {
        switch (quadPrec) {
            case Precision::VERYLOW: return 1;
            case Precision::LOW:     return 3;
            case Precision::MEDIUM:  return 7;
            case Precision::HIGH:    return 13;
        };
        } ();
}

/*
void Triangle::buildQuadCoeffs(Precision quadPrec) {

    switch (quadPrec) {
        case Precision::VERYLOW:
            quadCoeffs.reserve(1);

    }
}*/

void Triangle::buildQuads(Precision quadPrec) {

    auto baryPt = [&](double w0, double w1, double w2) {
        return w0*Xs[0] + w1*Xs[1] + w2*Xs[2];
        };

    //auto baryPt = [&](const vec3d& ws) {
    //    return ws[0]*Xs[0] + ws[1]*Xs[1] + ws[2]*Xs[2];
    //};

    // TODO: use permutation function
    switch (quadPrec) {
        case Precision::VERYLOW:
            quads.reserve(1);
            quads.emplace_back(center, 1.0/2.0);
            break;

        case Precision::LOW: {
            quads.reserve(3);

            constexpr double weight = 1.0/6.0;
            quads.emplace_back(baryPt(2.0/3.0, 1.0/6.0, 1.0/6.0), weight);
            quads.emplace_back(baryPt(1.0/6.0, 2.0/3.0, 1.0/6.0), weight);
            quads.emplace_back(baryPt(1.0/6.0, 1.0/6.0, 2.0/3.0), weight);
            break;
        }

        case Precision::MEDIUM: {
            quads.reserve(7);

            constexpr double weight0 = 0.1125;
            quads.emplace_back(center, weight0);

            constexpr double weight1 = 0.066197076394253;
            constexpr double alpha = 0.059715871789770, beta = 0.470142064105115;
            quads.emplace_back(baryPt(alpha, beta, beta), weight1);
            quads.emplace_back(baryPt(beta, alpha, beta), weight1);
            quads.emplace_back(baryPt(beta, beta, alpha), weight1);

            constexpr double weight2 = 0.0629695902724135;
            constexpr double gamma = 0.797426985353087, delta = 0.101286507323456;
            quads.emplace_back(baryPt(gamma, delta, delta), weight2);
            quads.emplace_back(baryPt(delta, gamma, delta), weight2);
            quads.emplace_back(baryPt(delta, delta, gamma), weight2);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);

            assert(quads.size() == 7);
            break;
        }

        case Precision::HIGH: {
            quads.reserve(13);

            constexpr double weight0 = -0.074785022233841;
            quads.emplace_back(center, weight0);

            constexpr double weight1 = 0.087807628716604;
            constexpr double alpha = 0.479308067841920, beta = 0.260345966079040;
            quads.emplace_back(baryPt(alpha, beta, beta), weight1);
            quads.emplace_back(baryPt(beta, alpha, beta), weight1);
            quads.emplace_back(baryPt(beta, beta, alpha), weight1);

            constexpr double weight2 = 0.026673617804419;
            constexpr double gamma = 0.869739794195568, delta = 0.065130102902216;
            quads.emplace_back(baryPt(gamma, delta, delta), weight2);
            quads.emplace_back(baryPt(delta, gamma, delta), weight2);
            quads.emplace_back(baryPt(delta, delta, gamma), weight2);

            constexpr double weight3 = 0.0385568804451285;
            constexpr double eps = 0.048690315425316, zeta = 0.312865496004874, theta = 0.638444188569810;
            quads.emplace_back(baryPt(eps, zeta, theta), weight3);
            quads.emplace_back(baryPt(theta, eps, zeta), weight3);
            quads.emplace_back(baryPt(zeta, theta, eps), weight3);
            quads.emplace_back(baryPt(eps, theta, zeta), weight3);
            quads.emplace_back(baryPt(zeta, eps, theta), weight3);
            quads.emplace_back(baryPt(theta, zeta, eps), weight3);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) + 6.0*weight3 - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);

            break;
        }
    }
}

/*
bool Triangle::isAdjacent(const std::shared_ptr<Triangle>& tri) {

}
*/