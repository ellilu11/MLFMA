#pragma once

#include <chrono>

using Time = std::chrono::duration<double, std::milli>;
using Clock = std::chrono::high_resolution_clock;

struct ClockTimes {
    ClockTimes() = default;

    Time M2X{ 0 };
    Time X2X{ 0 };
    Time X2L{ 0 };
    Time P2L{ 0 };
    Time L2L{ 0 };
    Time L2P{ 0 };
    Time M2P{ 0 };
    Time P2P{ 0 };
};

