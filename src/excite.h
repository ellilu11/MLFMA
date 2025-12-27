#pragma once

#include "types.h"

namespace Excitation {
    struct PlaneWave;

    // struct HertzDipole;
}

struct Excitation::PlaneWave {
    PlaneWave()
        : pol(vec3d{ 1,0,0 }), 
          wavevec(vec3d{ 0,0,1 }), 
          wavenum(wavevec.norm()), 
          amplitude(1.0)
    {};

    PlaneWave(
        const vec3d& pol,
        const vec3d& wavehat,
        double wavenum, double amplitude)
        : pol(pol), 
          wavevec(wavenum*wavehat), 
          wavenum(wavenum), 
          amplitude(amplitude)
    {};

    vec3d pol;          // unit polarization
    vec3d wavevec;      // wavevector
    double wavenum;     // wavenumber
    double amplitude;   // amplitude
};