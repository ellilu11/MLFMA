#pragma once

#include "interp.h"

struct Tables {
    Tables() = default;
    Tables(const int maxLevel,
        const int order,
        std::vector<realVec>& thetas, 
        std::vector<realVec>& phis,
        std::vector<int>& Ls,
        const double wavenum,
        const double rootLeng)
    {
        buildAngularTables(maxLevel, thetas, phis, wavenum);
        buildInterpThetaTable(maxLevel, order, thetas);
        buildInterpPhiTable(maxLevel, order, phis);
        buildTranslationTable(maxLevel, order, Ls, wavenum, rootLeng);

    }

    void buildAngularTables(
        const int, const std::vector<realVec>&, const std::vector<realVec>&, const double);

    void buildInterpThetaTable(const int, const int, const std::vector<realVec>&);

    void buildInterpPhiTable(const int, const int, const std::vector<realVec>&);

    void buildTranslationTable(
        const int, const int, const std::vector<int>&, const double, const double);

    // Angular tables
    std::vector<std::vector<mat3d>> ImKK;
    std::vector<std::vector<vec3d>> kvec;
    std::vector<std::vector<mat23d>> matToThPh;
    // std::vector<std::vector<mat3d>> matToThPh;

    // Lagrange interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<size_t>> ts;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<size_t>> ss;

    std::vector<std::vector<realVec>> interpPsi;

    // M2L translation tables
    std::vector<std::vector<cmplxVec>> transl;
    realVec normedDists;

};
