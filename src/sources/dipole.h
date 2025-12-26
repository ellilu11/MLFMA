#pragma once

#include <random>
#include "source.h"

class Dipole final : public Source {
public:
    Dipole() = default;

    Dipole(std::shared_ptr<PlaneWave>, const vec3d&);

    Dipole(std::shared_ptr<PlaneWave>, const vec3d&, const vec3d&);

    void buildVoltage() override;

    void buildCurrent() override;

    vec3d getCenter() const override { return pos; }

    friend std::ostream& operator<<(std::ostream& os, Dipole& src) {
        os << src.pos << ' ' << src.pol << '\n';
        return os;
    }

    /* getRadAlongDir(X,kvec)
     * Return the outgoing radiated amplitude at X along direction kvec
     * due to this dipole
     * X    : observation point (Cartesian)
     * kvec : wavevector
     */
    vec3cd getRadAlongDir(const vec3d& X, const vec3d& kvec) const override {
        return exp(iu*kvec.dot(X-pos)) * phat;
    }

    /* getRadAtPoint(X)
     * Return the radiated field due to this dipole
     * at field point X
     */
    //vec3cd getRadAtPoint(const vec3d& X) const override {
    //    assert(X != pos);
    //    return Math::dyadicG(X - pos, Einc->wavenum) * phat;
    //}

    /* getIntegratedRad(src)
     * Return the radiated field due to src tested with this dipole
     */
    cmplx getIntegratedRad(const std::shared_ptr<Source> src) const override {
        const auto srcDip = dynamic_pointer_cast<Dipole>(src);

        assert(pos != srcDip->pos);

        const auto& rad = Math::dyadicG(pos - srcDip->pos, Einc->wavenum) * srcDip->phat;

        return conj(rad.dot(phat));
    }

private:
    vec3d pos;  // position

    double pmag; // pol. density magnitude 
    vec3d pol; // pol. density vector
    vec3d phat; // unit pol. density vector

};