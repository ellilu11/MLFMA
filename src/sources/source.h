#pragma once

#include "../excite.h"
#include "../math.h"

class Source;

using SrcVec = std::vector<std::shared_ptr<Source>>;

class Source {

public:
    Source() = default;
    
    Source(std::shared_ptr<PlaneWave> Einc)
        : Einc(std::move(Einc)), rhs(0.0), current(0.0), sol(0.0) 
    {};

    cmplx getCurrent() const { return current; }

    cmplx getSol() const { return sol; }

    void addToSol(cmplx sol_) { sol += sol_; }

    void resetSol() { sol = 0.0; }

    void printSol(std::ofstream& f) const {
        f << sol.real() << ' ' << sol.imag() << '\n';
    }

    virtual vec3d getCenter() const = 0;

    virtual void buildRHS() = 0;

    virtual void buildCurrent() = 0;

    virtual vec3cd getRadAlongDir(const vec3d&, const vec3d&) const = 0;

    virtual vec3cd getIncAlongDir(const vec3d&, const vec3d&) const = 0;

    virtual vec3cd getRadAtPoint(const vec3d&) const = 0;

    virtual cmplx getIntegratedRad(const std::shared_ptr<Source>) const = 0;

    friend std::ostream& operator<<(std::ostream& os, Source& src) {
        os << src.getCenter() << '\n';
        return os;
    }

protected:
    std::shared_ptr<PlaneWave> Einc;
    cmplx rhs;
    cmplx current;
    cmplx sol;

};