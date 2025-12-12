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
    Tables(bool isDefault)
    {
        buildAngularTables();
        
        buildInterpThetaTable();
        buildInterpPhiTable();
        
        buildTranslationTable();
        // buildInterpPsiTable();
    }
    
    void buildAngularTables();

    void buildInterpThetaTable();

    void buildInterpPhiTable();

    void buildTranslationTable();

    void buildInterpPsiTable();

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

