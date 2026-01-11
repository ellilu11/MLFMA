#pragma once

#include "math.h"

namespace Interp {
    std::pair<realVec, realVec> gaussLegendre(
        const int, const double EPS, const double, const double);

    int getNearGLNodeIdx(
        const double, const int, const double, const double);

    double evalLagrangeBasis(const double, const realVec&, const int);

    double evalTrigBasis(const double, const realVec&, const int);
}