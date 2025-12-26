#include "dipole.h"

Dipole::Dipole(
    std::shared_ptr<PlaneWave> Einc, const vec3d& X)
    : Source(Einc), pos(X), pmag(Phys::p0), pol(vec3d(pmag, 0, 0)), phat(pol/pmag)
{
    // buildVoltage();

    buildCurrent();
};

Dipole::Dipole(
    std::shared_ptr<PlaneWave> Einc, const vec3d& X, const vec3d& P)
    : Dipole(Einc, X)
{
    pol = P;
    pmag = P.norm();
    phat = P / pmag;

    // std::cout << "Pol: " << P << '\n';
};

void Dipole::buildVoltage() {

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    voltage = -Einc->amplitude * exp(iu*kvec.dot(pos)) 
                * pol.dot(Einc->pol);

}

void Dipole::buildCurrent() {

    current = iu * Phys::c0 * Einc->wavenum * pmag; // |J| = i \omega |P|

    // std::cout << current << '\n';
}







