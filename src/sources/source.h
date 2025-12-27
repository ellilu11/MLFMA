#pragma once

#include "../excite.h"
#include "../math.h"

class Source;

using SrcVec = std::vector<std::shared_ptr<Source>>;

class Source {

public:
    Source() = default;
    
    Source(std::shared_ptr<Excitation::PlaneWave> Einc, size_t idx)
        : Einc(std::move(Einc)), idx(idx), voltage(0.0)
    {}

    cmplx getVoltage() const { return voltage; }

    //cmplx getCurrent() const { return current; }

    //cmplx getSol() const { return sol; }

    //void addToSol(cmplx sol_) { sol += sol_; }

    //void resetSol() { sol = 0.0; }

    //void printSol(std::ofstream& f) const {
    //    f << sol.real() << ' ' << sol.imag() << '\n';
    //}

    size_t getIdx() const { return idx; }

    virtual vec3d getCenter() const = 0;

    virtual void buildVoltage() = 0;

    // virtual void buildCurrent() = 0;

    virtual vec3cd getRadAlongDir(const vec3d&, const vec3d&) const = 0;

    virtual cmplx getIntegratedRad(const std::shared_ptr<Source>) const = 0;

protected:
    std::shared_ptr<Excitation::PlaneWave> Einc;
    cmplx voltage;
    //cmplx current;
    //cmplx sol;

    size_t idx;

};