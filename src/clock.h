#pragma once

#include <chrono>

using Time = std::chrono::duration<double, std::milli>;
using Clock = std::chrono::high_resolution_clock;

struct ClockTimes {
    ClockTimes() = default;

    Time S2M{ 0 };
    Time M2M{ 0 };
    Time M2L{ 0 };
    Time L2L{ 0 };
    Time L2T{ 0 };
    Time S2T{ 0 };

    Time M2L_lookup{ 0 };
};