#pragma once

#include "math.h"

namespace Phys {
    constexpr double c0 = 299792458.0;
    constexpr double mu0 = 1.256637E-6;
    constexpr double p0 = 1.0E-8; // Dipole moment

    const cmplx C = -iu * c0 * mu0 / (4.0 * PI);
}