#pragma once

#include "clock.h"
#include "interp.h"
#include "maps.h"

class Tables {

public: 
    Tables() = default;

    Tables(int maxLevel, int order)
        : maxLevel(maxLevel),
          order(order),
          dists(Math::getINodeDistances()),
          rhats(Math::getINodeDirections())
    {
        buildAngularTables();
        
        buildInterpTables();
        
        buildTranslationTable();
    }
    
    // Angular tables
    std::vector<std::vector<vec3d>> khat;
    std::vector<std::vector<mat23d>> toThPh;
    std::vector<std::vector<mat3d>> ImRR;

    // M2M interpolation tables
    std::vector<std::vector<interpPair>> interpTheta;
    std::vector<std::vector<interpPair>> interpPhi;

    // L2L interpolation tables
    std::vector<std::vector<interpPair>> invInterpTheta;
    std::vector<std::vector<interpPair>> invInterpPhi;

    // M2L translation table
    std::vector<VecHashMap<vecXcd>> transl;

private:
    std::vector<interpPair> getInterpThetaAtLvl(int, int);

    std::vector<interpPair> getInterpPhiAtLvl(int, int);

    Map<vecXcd> getAlphaAtLvl(int);

    HashMap<interpPair> getInterpPsiAtLvl(int);

    void buildAngularTables();

    void buildInterpTables();

    void buildTranslationTable();

    int maxLevel;
    int order;
    realVec dists;
    std::vector<vec3d> rhats;
};

