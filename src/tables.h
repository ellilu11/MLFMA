#pragma once

#include "math.h"

struct Tables {
    Tables() = default;
    Tables(const int maxLevel, const double wavenum, std::vector<realVec>& thetas) {
        buildDirTables(maxLevel, wavenum, thetas);
        // buildYlmTables(order);
    }

    void buildDirTables(const int, const double, const std::vector<realVec>&);
    void buildInterpThetaTable(const int, const int, const std::vector<realVec>&);
    void buildInterpPhiTable(const int, const int, const std::vector<realVec>&);
    // void buildYlmTables(const int);

    // Directional tables
    std::vector<std::vector<mat3d>> ImKK;
    std::vector<std::vector<vec3d>> kvec;

    // Lagrange interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<size_t>> T;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<size_t>> S;

    // Ylm tables
    /*std::vector<realVec> coeffYlm_;
    std::vector<realVec> fallingFact_;
    std::vector<realVec> legendreSum_;
    std::vector<realVec> fracCoeffYlm_;
    std::vector<realVec> A_;
    std::vector<realVec> Aexp_;*/
};
