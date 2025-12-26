#pragma once

#include <map>
#include <unordered_map>

namespace Maps {

    constexpr double DistEPS = 1.0E-3;
    constexpr double PsiEPS = 1.0E-3; // Pick largest value that avoids collisions

    struct Comp {
        bool operator()(double x, double y) const noexcept {
            return x + DistEPS < y;
        }
    };

    struct HashFunc {
        std::size_t operator()(double x) const noexcept {
            long long q = static_cast<long long>(std::llround(x / PsiEPS));
            return std::hash<long long>{}(q);
        }
    };

    struct HashEq {
        bool operator()(double x, double y) const noexcept {
            return std::fabs(x - y) <= PsiEPS;
        }
    };
}

template <typename T, typename U>
using Map = std::map<T, U, Maps::Comp>;

template <typename T, typename U>
using HashMap = std::unordered_map<T, U, Maps::HashFunc, Maps::HashEq>;