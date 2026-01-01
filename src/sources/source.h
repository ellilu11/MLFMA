#pragma once

#include "../excite.h"
#include "../math.h"

class Source;

using SrcVec = std::vector<std::shared_ptr<Source>>;

class Source {

public:
    Source() = default;
    
    Source(std::shared_ptr<Excitation::PlaneWave> Einc)
        : Einc(std::move(Einc)), voltage(0.0), current(0.0), sol(0.0) 
    {};

    cmplx getVoltage() const { return voltage; }

    cmplx getSol() const { return sol; }

    void setIdx(size_t idx_) { idx = idx_; }

    virtual vec3d getCenter() const = 0;

    virtual void buildVoltage() = 0;

    virtual vec3cd getRadAlongDir(const vec3d&, const vec3d&) const = 0;

    virtual cmplx getIntegratedRad(const std::shared_ptr<Source>) const = 0;

protected:
    std::shared_ptr<Excitation::PlaneWave> Einc;
    cmplx voltage;
    size_t idx;

};