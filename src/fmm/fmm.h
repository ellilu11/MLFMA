#pragma once

#include "clock.h"
#include "types.h"

extern ClockTimes t;

namespace FMM {
    constexpr int DIM = 3;
    constexpr int numDir = 26;

    constexpr double EPS_NR = 1.0E-9; // Newton-Raphson precision

    enum class Dir {
        W, E, S, N, D, U,
        SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
        DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
    };

    struct Angles;

    class Node;
    class Stem;
    class Leaf;

    class Tables;
}