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
    //for (const auto& X : Xs)
    //    std::cout << X << '\n';
    //std::cout << '\n';

    for (int i = 0; i < 3; ++i) {
        const int ipp = Math::wrapIdxToRange(i+1, 3);
        Ds[i] = Xs[ipp] - Xs[i];

        // const int imm = Math::wrapIdxToRange(i-1, 3);
        // Dms[i] = Xs[i] - Xs[imm];
    }

    buildQuads(quadPrec);
};

int Triangle::prec2Int(const Precision quadPrec) {
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

void Triangle::buildQuads(const Precision quadPrec) {

    auto baryPt = [&](double w0, double w1, double w2) {
        return w0*Xs[0] + w1*Xs[1] + w2*Xs[2];
    };

    switch (quadPrec) {
        case Precision::VERYLOW:
            quads.emplace_back(center, 1.0/2.0);
            break;

        case Precision::LOW: {
            constexpr double weight = 1.0/6.0;
            quads.emplace_back(baryPt(2.0/3.0, 1.0/6.0, 1.0/6.0), weight);
            quads.emplace_back(baryPt(1.0/6.0, 2.0/3.0, 1.0/6.0), weight);
            quads.emplace_back(baryPt(1.0/6.0, 1.0/6.0, 2.0/3.0), weight);
            break;
        }

        case Precision::MEDIUM: {
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

            break;
        }

        case Precision::HIGH: {
            // TODO: 13-point quadrature
        }
    }
}

bool Triangle::isAdjacent(const std::shared_ptr<Triangle>& tri) {

    return false;
}