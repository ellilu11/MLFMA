#pragma once

#include <map>
#include <unordered_map>

namespace Maps {

    // Pick largest values that avoid collisions
    constexpr double EPS_dr = 1.0E-1;
    constexpr double EPS_psi = 1.0E-9;
    constexpr double EPS_dX = 1.0E-3;

    struct DoubleComp {
        bool operator()(double x, double y) const noexcept {
            return x + EPS_dr < y;
        }
    };

    struct DoubleFunc {
        std::size_t operator()(double x) const noexcept {
            long long q = static_cast<long long>(std::llround(x / EPS_psi));
            return std::hash<long long>{}(q);
        }
    };

    struct DoubleEq {
        bool operator()(double x, double y) const noexcept {
            return std::fabs(x - y) <= EPS_psi;
        }
    };

    struct VecFunc {
        std::size_t operator()(const vec3d& X) const noexcept {
            /*auto q = [](double x) {
                return static_cast<long long>(std::llround(x / EPS));
            };
            size_t h = 0;
            hash_combine(h, q(X[0]));
            hash_combine(h, q(X[1]));
            hash_combine(h, q(X[2]));
            return h;*/

            // Puzzle: Why does hashing the norm work?
            long long q = static_cast<long long>(std::llround(X.norm() / EPS_dX));
            return std::hash<long long>{}(q);
        }
    };

    struct VecEq {
        bool operator()(const vec3d& X, const vec3d& Y) const noexcept {
            return ((X-Y).norm()) <= EPS_dX;
        }
    };
}

// TODO: Concepts
template <typename T>
using Map = std::map<double, T, Maps::DoubleComp>;

template <typename T>
using HashMap = std::unordered_map<double, T, Maps::DoubleFunc, Maps::DoubleEq>;

template <typename T>
using VecHashMap = std::unordered_map<vec3d, T, Maps::VecFunc, Maps::VecEq>;