#include "triangle.h"

using namespace std;

Triangle::Triangle(
    const vec3i& vIdx,
    const std::vector<vec3d>& vList,
    const Precision quadPrec)
    : vIdx(vIdx),
    Xs({ vList[vIdx[0]], vList[vIdx[1]], vList[vIdx[2]] })
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
            case Precision::LOW:     return 1;
            case Precision::MEDIUM:  return 3;
            case Precision::HIGH:    return 7;
        };
        } ();
}

void Triangle::buildQuads(const Precision quadPrec) {

    switch (quadPrec) {
        case Precision::LOW :
            quadNodes.push_back((Xs[0] + Xs[1] + Xs[2]) / 3.0);
            quadWeight = 1.0/2.0;
            break;

        case Precision::MEDIUM :
            quadNodes.resize(3);
            quadNodes[0] = 2.0/3.0*Xs[0] + 1.0/6.0*Xs[1] + 1.0/6.0*Xs[2];
            quadNodes[1] = 1.0/6.0*Xs[0] + 2.0/3.0*Xs[1] + 1.0/6.0*Xs[2];
            quadNodes[2] = 1.0/6.0*Xs[0] + 1.0/6.0*Xs[1] + 2.0/3.0*Xs[2];
            quadWeight = 1.0/6.0;
            break;

        case Precision::HIGH :
            quadNodes.resize(7);
            // TODO: 7-point quadrature
            break;
    }
}

bool Triangle::isAdjacent(const std::shared_ptr<Triangle>& tri) {



    return false;
}