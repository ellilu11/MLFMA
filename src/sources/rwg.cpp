#include "rwg.h"

RWG::RWG(
    std::shared_ptr<PlaneWave> Einc,
    const Eigen::Vector4i& idx,
    const std::vector<vec3d>& vertices,
    const TriVec& triangles)
    : Source(Einc),
      v0(vertices[idx[0]]), v1(vertices[idx[1]]),
      triPlus(triangles[idx[2]]), triMinus(triangles[idx[3]]),
      center((v0+v1)/2.0)

{
    for (const auto& vIdx : triPlus->getVidx())
        if (vIdx != idx[0] && vIdx != idx[1])
            vPlus = vertices[vIdx];

    for (const auto& vIdx : triMinus->getVidx())
        if (vIdx != idx[0] && vIdx != idx[1])
            vMinus = vertices[vIdx];

    leng = (v0-v1).norm();

    // buildRHS();

    buildCurrent();

    // std::cout << '(' << v0 << ") (" << v1 << ") ("
    //    << vPlus << ") (" << vMinus << ") " << leng << '\n';
};

void RWG::buildRHS() {

    cmplx rhs(0,0);

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        rhs += weightPlus * exp(iu*kvec.dot(quadNode)) 
                    * (quadNode - vPlus).dot(Einc->pol);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        rhs += weightMinus * exp(iu*kvec.dot(quadNode))
                    * (vMinus - quadNode).dot(Einc->pol);

    rhs *= -Einc->amplitude * leng;
    
}

void RWG::buildCurrent() {
    current = 1.0;

    // TODO: Predict current from rhs vector

}

/* getRadAlongDir(X,kvec)
 * Return the outgoing radiated amplitude due to this
 * RWG at X along direction kvec
 * X    : observation point (Cartesian)
 * kvec : wavevector 
 */ 
vec3cd RWG::getRadAlongDir(
    const vec3d& X, const vec3d& kvec) const {
    
    vec3cd rad = vec3cd::Zero();

    auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        rad += weightPlus * exp(iu*kvec.dot(X-quadNode))
                * (vPlus - quadNode);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        rad += weightMinus * exp(iu*kvec.dot(X-quadNode))
                * (quadNode - vMinus);

    return current * leng * rad;
}

/* getIncAlongDir(X,kvec)
 * Return the incoming radiated amplitude at this
 * RWG due to field at X along direction kvec
 * X    : source point (Cartesian)
 * kvec : wavevector
 */
vec3cd RWG::getIncAlongDir(
    const vec3d& X, const vec3d& kvec) const {

    vec3cd inc = vec3cd::Zero();

    auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        inc += weightPlus * exp(iu*kvec.dot(quadNode-X))
                * (vPlus - quadNode);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        inc += weightMinus * exp(iu*kvec.dot(quadNode-X))
                * (quadNode - vMinus);

    return leng * inc; // no current term here!
}

/* getRadAtPoint(X)
 * Return the full radiated field with given wavenum
 * due to this RWG at X
 */
vec3cd RWG::getRadAtPoint(const vec3d& X) const {

    vec3cd rad = vec3cd::Zero();

    const double wavenum = Einc->wavenum;
   
    auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus) {

        const auto& dyadic = Math::dyadicG(X - quadNode, wavenum);

        rad += weightPlus * dyadic * (vPlus - quadNode);
    }

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus) {
    
        const auto& dyadic = Math::dyadicG(X - quadNode, wavenum);
        
        rad += weightMinus * dyadic * (quadNode - vMinus);
    }

    return current * leng * rad;
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx RWG::getIntegratedRad(const std::shared_ptr<Source> src) const {

    cmplx integratedRad = 0.0;

    const double wavenum = Einc->wavenum;

    auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        integratedRad += weightPlus * src->getRadAtPoint(quadNode).dot(vPlus - quadNode);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        integratedRad += weightMinus * src->getRadAtPoint(quadNode).dot(quadNode - vMinus);

    return leng * integratedRad;
    
}