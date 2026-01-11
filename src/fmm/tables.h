#pragma once

#include "../interp.h"
#include "../maps.h"
#include "fmm.h"

class FMM::Tables {

public:
    Tables() = default;

    Tables(int level, int maxLevel)
        : level(level)
    {
        buildAngularTables();
        if (level < maxLevel) buildInterpTables();
        buildTranslationTable();
    }

    static void buildDists() {
        dists = Math::getINodeDistances();
        rhats = Math::getINodeDirections();
        dXs = Math::getINodeDistVecs();
    }

    // Angular tables
    std::vector<vec3d> khat;
    std::vector<mat23d> toThPh;
    std::vector<mat3d> ImRR;

    // M2M interpolation tables
    std::vector<interpPair> interpTheta;
    std::vector<interpPair> interpPhi;

    // L2L interpolation tables
    std::vector<interpPair> invInterpTheta;
    std::vector<interpPair> invInterpPhi;

    // M2L translation table
    VecHashMap<vecXcd> transl;

private:
    std::vector<interpPair> getInterpTheta(int, int);

    std::vector<interpPair> getInterpPhi(int, int);

    Map<vecXcd> getAlpha();

    HashMap<interpPair> getInterpPsi();

    void buildAngularTables();

    void buildInterpTables();

    void buildTranslationTable();

    static realVec dists;
    static std::vector<vec3d> rhats;
    static std::array<vec3d,316> dXs;

    int level;
};
