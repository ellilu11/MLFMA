#pragma once

#include <random>
#include "source.h"

class Dipole final : public Source {
public:
    Dipole() = default;

    Dipole(std::shared_ptr<PlaneWave>, const vec3d&);

    void buildVoltage() override;

    void buildCurrent() override;

    vec3d getCenter() const override { return pos; }

    /* getRadAlongDir(X,kvec)
     * Return the outgoing radiated amplitude at X along direction kvec
     * due to this dipole
     * X    : observation point (Cartesian)
     * kvec : wavevector
     */
    vec3cd getRadAlongDir(const vec3d& X, const vec3d& kvec) const override {
        return exp(iu*kvec.dot(X-pos)) * phat;
    }

    /* getIncAlongDir(X,kvec)
     * Return the incoming radiated amplitude from X along direction kvec
     * at this dipole
     * X    : source point (Cartesian)
     * kvec : wavevector
     */
    //vec3cd getIncAlongDir(const vec3d& X, const vec3d& kvec) const override {
    //    return conj(exp(iu*kvec.dot(pos-X))) * phat;
    //}

    /* getRadAtPoint(X)
     * Return the radiated field due to this dipole
     * at field point X
     */
    vec3cd getRadAtPoint(const vec3d& X) const override {
        assert(X != pos);
        return Math::dyadicG(X - pos, Einc->wavenum) * phat;
    }

    /* getIntegratedRad(src)
     * Return the radiated field due to src tested with this dipole
     */
    cmplx getIntegratedRad(const std::shared_ptr<Source> src) const override {
        const auto srcDip = dynamic_pointer_cast<Dipole>(src);

        return conj(srcDip->getRadAtPoint(pos).dot(phat));
    }

private:
    vec3d pos;  // position
    vec3d phat; // unit pol. density vector
    double pol; // pol. density magnitude 
};