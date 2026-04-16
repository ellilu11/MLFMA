#pragma once

#include "triangle.h"
#include "../types/types.h"

namespace Mesh {
    struct Mesh::TriPairs {

    public:
        TriPairs() = default;

        TriPairs(size_t);

        void clear() { *this = {}; }

        static std::pair<Triangle, Triangle> getTriPair(pair2i iTris) {
            return { glTris[iTris.first], glTris[iTris.second] };
        }

    private:
        void buildNumCommon();

        void buildMomentsEFIE();

        void buildMomentsMFIE();

        void buildIntegratedInvR();

        void buildIntegratedInvRcubed();

    public:
        std::vector<MomentsEFIE> momentsEFIE;
        std::vector<MomentsMFIE> momentsMFIE;
        std::vector<MomentsMFIE> momentsMFIE2; // symmetric case

        std::vector<intRads> intsInvR;
        std::vector<intRads> intsInvR2; // symmetric case

        std::vector<intRads> intsInvRcubed;
        std::vector<intRads> intsInvRcubed2; // symmetric case

        std::vector<int> nCommons;  // number of common vertices

    private:
        std::vector<std::pair<pair2i, size_t>> pairsToIdx;
        size_t nPair; // number of triangle pairs
    };

    TriPairs glTriPairs; // triangle pairs data
}

