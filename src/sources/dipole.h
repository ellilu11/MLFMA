#pragma once

#include <random>
#include "source.h"

class Dipole final : public Source {
public:
    Dipole() = default;

    Dipole(
        std::shared_ptr<Excitation::PlaneWave> Einc, size_t idx, const vec3d& X)
        : Source(Einc, idx), pos(X), 
        pmag(Phys::p0), pol(vec3d(pmag, 0, 0)), phat(pol/pmag)
    {
        buildVoltage();
    };

    Dipole(
        std::shared_ptr<Excitation::PlaneWave> Einc, size_t idx, const vec3d& X, const vec3d& P)
        : Dipole(Einc, idx, X)
    {
        pol = P;
        pmag = P.norm();
        phat = P / pmag;

        // std::cout << "Pol: " << P << '\n';
    };

    void buildVoltage() override {
        voltage = -Einc->amplitude * exp(iu*Einc->wavevec.dot(pos))
            * pol.dot(Einc->pol);
    }

    //void buildCurrent() override {
    //    current = iu * Phys::c0 * Einc->wavenum * pmag; // |J| = i \omega |P|
    //    // std::cout << current << '\n';
    //}

    vec3d getCenter() const override { return pos; }

    /* getRadAlongDir(X,kvec)
     * Return the outgoing radiated amplitude at X along direction kvec
     * due to this dipole
     * X    : observation point (Cartesian)
     * kvec : wavevector
     */
    vec3cd getRadAlongDir(
        const vec3d& X, const vec3d& kvec) const override 
    {
        return exp(iu*kvec.dot(X-pos)) * phat;
    }

    /* getIntegratedRad(src)
     * Return the radiated field due to src tested with this dipole
     */
    cmplx getIntegratedRad(const std::shared_ptr<Source> src) const override {
        const auto srcDip = dynamic_pointer_cast<Dipole>(src);

        if (pos == srcDip->pos) return 0.0; // TODO: Radiation reaction field

        const auto& rad = Math::dyadicG(pos - srcDip->pos, Einc->wavenum) * srcDip->phat;

        return conj(rad.dot(phat));
    }

    friend std::ostream& operator<<(std::ostream& os, Dipole& src) {
        os << src.pos << ' ' << src.pol << '\n';
        return os;
    }

private:
    vec3d pos;  // position

    double pmag; // pol. density magnitude 
    vec3d pol; // pol. density vector
    vec3d phat; // unit pol. density vector

};