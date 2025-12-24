#include "dipole.h"

Dipole::Dipole(
    std::shared_ptr<PlaneWave> Einc, 
    const vec3d& X) // TODO: Pass in pol. density vector
    : Source(Einc), pos(X), phat(vec3d(1, 0, 0)), pol(1.0E-8)
{
    // buildVoltage();

    buildCurrent();
};

void Dipole::buildVoltage() {

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    voltage = -Einc->amplitude * exp(iu*kvec.dot(pos)) 
                * pol * phat.dot(Einc->pol);

}

void Dipole::buildCurrent() {

    current = iu * c0 * Einc->wavenum * pol; // |J| = i \omega |P|

    // std::cout << current << '\n';
}







