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

    size_t getIdx() const { return idx; }

    virtual vec3d getCenter() const = 0;

    virtual void buildVoltage() = 0;

    virtual vec3cd getRadAlongDir(const vec3d&, const vec3d&) const = 0;

    virtual cmplx getIntegratedRad(const std::shared_ptr<Source>) const = 0;

protected:
    std::shared_ptr<Excitation::PlaneWave> Einc;
    cmplx voltage;
    size_t idx;

};