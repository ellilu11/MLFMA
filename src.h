#pragma once

#include "math.h"

struct Src {
    Src()
        : amplitude(1.0),
        pol(vec3d{ 1,0,0 }),
        wavevec(vec3d{ 0,0,1 })
    {
    };

    double amplitude;   // amplitude
    vec3d pol;          // polarization (normalized)
    vec3d wavevec;      // wavevector
};