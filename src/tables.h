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

        // buildInterpPsiTable();
    }
    
    void buildAngularTables();

    std::vector<interpPair> getInterpThetaAtLvl(int, int);

    std::vector<interpPair> getInterpPhiAtLvl(int, int);

    void buildInterpTables();

    void buildTranslationTable();

    void buildInterpPsiTable();

    int order;

    // Angular tables
    std::vector<std::vector<vec3d>> khat;
    std::vector<std::vector<mat23d>> toThPh;

    // M2M interpolation tables
    std::vector<std::vector<interpPair>> interpTheta;
    std::vector<std::vector<interpPair>> interpPhi;

    // L2L interpolation tables
    std::vector<std::vector<interpPair>> invInterpTheta;
    std::vector<std::vector<interpPair>> invInterpPhi;

    // M2L translation tables
    std::vector<Map<double,vecXcd>> transl;
    std::vector<HashMap<double,interpPair>> interpPsi;

};

