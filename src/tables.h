#pragma once

#include "clock.h"
#include "interp.h"
#include "map.h"

struct Tables {
    Tables() = default;
    Tables(const Config& config)
        : order(config.interpOrder)
    {
        buildAngularTables();
        
        buildInterpTables();
        
        buildTranslationTable();
        buildInterpPsiTable();
    }
    
    void buildAngularTables();

    std::pair<std::vector<realVec>, std::vector<int>> 
        getInterpThetaAtLvl(int, int);

    std::pair<std::vector<realVec>, std::vector<int>> 
        getInterpPhiAtLvl(int, int);

    void buildInterpTables();

    void buildTranslationTable();
    void buildInterpPsiTable();

    int order;

    // Angular tables
    std::vector<std::vector<mat3d>> ImKK;
    std::vector<std::vector<vec3d>> khat;

    std::vector<std::vector<mat23d>> matToThPh;
    std::vector<std::vector<mat32d>> matFromThPh;

    // M2M interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<int>> idxTheta;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<int>> idxPhi;

    // L2L interpolation tables
    std::vector<std::vector<realVec>> invInterpTheta;
    std::vector<std::vector<int>> invIdxTheta;

    std::vector<std::vector<realVec>> invInterpPhi;
    std::vector<std::vector<int>> invIdxPhi;

    // M2L translation tables
    std::vector<Map<double,vecXcd>> transl;
    std::vector<HashMap<double,vecXcd>> interpPsi;
    std::vector<HashMap<double,int>> idxPsi;

};

