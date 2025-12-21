#include "dipole.h"

Dipole::Dipole(
    std::shared_ptr<PlaneWave> Einc, 
    const vec3d& X) // TODO: Pass in pol. density vector
    : Source(Einc), pos(X), phat(vec3d(1, 0, 0)), pol(1.0E-8)
{
    // buildRHS();

    buildCurrent();
};

void Dipole::buildRHS() {

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    rhs = -Einc->amplitude * exp(iu*kvec.dot(pos)) 
            * pol * phat.dot(Einc->pol);

}

void Dipole::buildCurrent() {

    current = iu * c0 * Einc->wavenum * pol; // |J| = i \omega |P|

    // std::cout << current << '\n';
}

/* getRadAlongDir(X,kvec)
 * Return the outgoing radiated amplitude at X along direction kvec 
 * due to this dipole 
 * X    : observation point (Cartesian)
 * kvec : wavevector
 */
vec3cd Dipole::getRadAlongDir(
    const vec3d& X, const vec3d& kvec) const {

    return current * exp(iu*kvec.dot(X-pos)) * phat; // TODO: Multiply in current later
}

/* getIncAlongDir(X,kvec)
 * Return the incoming radiated amplitude from X along direction kvec
 * at this dipole 
 * X    : source point (Cartesian)
 * kvec : wavevector
 */
vec3cd Dipole::getIncAlongDir(
    const vec3d& X, const vec3d& kvec) const {

    return conj(exp(iu*kvec.dot(pos-X))) * phat; 
}

/* getRadAtPoint(X)
 * Return the radiated field due to this dipole 
 * at field point X
 */
vec3cd Dipole::getRadAtPoint(const vec3d& X) const {
    assert(X != pos);

    const auto& dyadic = Math::dyadicG(X - pos, Einc->wavenum);

    return current * dyadic * phat; // TODO: Multiply in current later
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this dipole
 */
cmplx Dipole::getIntegratedRad(const std::shared_ptr<Source> src) const {

    const auto srcDip = dynamic_pointer_cast<Dipole>(src);

    return conj(srcDip->getRadAtPoint(pos).dot(phat));
}