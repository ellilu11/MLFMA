#pragma once

#include "../clock.h"
#include "../config.h"
#include "../exct.h"
#include "../math.h"
#include "../phys.h"
#include "../types/types.h"

extern const Config config;
extern auto t = ClockTimes();

namespace FMM {
    // Constants
    constexpr int DIM = 3;
    constexpr int numDir = 26;

    // Types
    class Node;
    using NodeVec = std::vector<std::shared_ptr<Node>>;
    using NodePair = std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>>;

    class Farfield;
    class Nearfield;

    struct Coeffs;
    struct Level;

    using interpPair = std::pair<std::vector<double>, int>;

    enum class Dir {
        W, E, S, N, D, U,
        SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
        DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
    };

    // Global data
    NodeVec glNodes;
    NodeVec glLeaves;
    std::vector<NodeVec> glNodesByLvl;
    std::vector<NodePair> glNonNearPairs;
    inline int maxLevel = 0;
    inline int glSrcIdx = 0;

    // Functions
    void reset() {
        glNodes.clear();
        glLeaves.clear(); 
        glNodesByLvl.clear();
        glNonNearPairs.clear();
        maxLevel = 0;
        glSrcIdx = 0;
    }
}