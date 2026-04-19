#pragma once

#include <chrono>

using Time = std::chrono::duration<double, std::milli>;
using Clock = std::chrono::high_resolution_clock;

struct ClockTimes {
    ClockTimes() = default;

    void printTimes() const {
        std::cout << "   S2M: " << S2M.count() << " ms\n"
                  << "   M2M: " << M2M.count() << " ms\n"
                  << "   M2L: " << M2L.count() << " ms\n"
                  << "   L2L: " << L2L.count() << " ms\n"
                  << "   L2T: " << L2T.count() << " ms\n"
                  << "   S2T: " << S2T.count() << " ms\n"
                  << "   ILU: " << ILU.count() << " ms\n";
    }

    void resetTimes() { *this = {}; }

    Time S2M{};
    Time M2M{};
    Time M2L{};
    Time L2L{};
    Time L2T{};
    Time S2T{};
    Time ILU{};
};