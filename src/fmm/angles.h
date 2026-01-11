#pragma once

#include "fmm.h"

struct FMM::Angles {
    Angles() = default;

    Angles(int level) {
        buildAngularSamples(level);

        buildAngularMatrices();
    }

    void buildAngularSamples(int);

    void buildAngularMatrices();

    void printAngles(std::ofstream&, std::ofstream&);

    pair2i getNumAngles() const {
        return std::make_pair(thetas.size(), phis.size());
    }

    size_t getNumAllAngles() const {
        return thetas.size() * phis.size();
    }

    std::vector<vec3d> khat;
    std::vector<mat23d> toThPh;
    std::vector<mat3d> ImRR;

    realVec thetas;  // theta samples
    realVec weights; // weights of theta samples
    realVec phis;    // phi samples
    int L;           // M2L series truncation number
};