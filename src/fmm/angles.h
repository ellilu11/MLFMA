#pragma once

#include "fmm.h"

struct FMM::Angles {
    Angles() = default;

    Angles(int);

    void printAngles(std::ofstream&, std::ofstream&);

    pair2i getNumAngles() const {
        return std::make_pair(thetas.size(), phis.size());
    }

    size_t getNumAllAngles() const {
        return thetas.size() * phis.size();
    }

    realVec thetas;  // theta samples
    realVec weights; // weights of theta samples
    realVec phis;    // phi samples
    int L;           // M2L series truncation number
};