#pragma once

#include <map>
#include "interp.h"

constexpr double EPS = 1.0E-2;

struct Comp {
    bool operator()(double x, double y) const {
        return x + EPS < y;
    }
};

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
        // buildInterpPsiTable(thetas, phis, maxLevel, order, wavenum);

    }
    
    void buildAngularTables(
        const int, const std::vector<realVec>&, const std::vector<realVec>&, const double);

    void buildInterpThetaTable(const int, const int, const std::vector<realVec>&);

    void buildInterpPhiTable(const int, const int, const std::vector<realVec>&);

    void buildTranslationTable(
        const int, const int, const std::vector<int>&, const double, const double);

    void buildInterpPsiTable(
        const std::vector<realVec>&, const std::vector<realVec>&, int, int, double);

    // Angular tables
    std::vector<std::vector<mat3d>> ImKK;
    std::vector<std::vector<vec3d>> kvec;
    std::vector<std::vector<Eigen::Matrix<double,2,3>>> matToThPh;
    std::vector<std::vector<Eigen::Matrix<double,3,2>>> matFromThPh;
    
    std::vector<std::vector<mat3d>> matToSph;
    std::vector<std::vector<mat3d>> matFromSph;

    // Lagrange interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<int>> ts;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<int>> ss;

    // M2L translation tables
    std::vector<std::map<double,cmplxVec,Comp>> transl;
    std::vector<std::map<double,realVec,Comp>> interpPsi;

};

