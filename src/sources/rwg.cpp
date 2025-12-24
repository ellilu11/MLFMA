#include "rwg.h"

RWG::RWG(
    std::shared_ptr<PlaneWave> Einc,
    const Eigen::Vector4i& idx,
    const std::vector<vec3d>& vertices,
    const TriVec& triangles)
    : Source(Einc),
      tris({triangles[idx[2]], triangles[idx[3]]}),
      X0(vertices[idx[0]]), 
      X1(vertices[idx[1]]),
      center((X0+X1)/2.0), 
      leng((X0-X1).norm())
{

    // Find non-common vertices
    for (int i = 0; i < 2; ++i)
        for (const auto& vIdx : tris[i]->vIdx)
            if (vIdx != idx[0] && vIdx != idx[1])
                Xpm[i] = vertices[vIdx];

    // buildVoltage();

    buildCurrent();

    //std::cout << '(' << X0 << ") (" << X1 << ") ("
    //    << Xpm[0] << ") (" << Xpm[1] << ") " << leng << '\n';
};

void RWG::buildVoltage() {

    cmplx voltage(0,0);

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    /*auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        rhs += weightPlus * exp(iu*kvec.dot(quadNode)) 
                    * (quadNode - vPlus).dot(Einc->pol);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        rhs += weightMinus * exp(iu*kvec.dot(quadNode))
                    * (vMinus - quadNode).dot(Einc->pol);
    */

    voltage *= -Einc->amplitude * leng;
    
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
   
    const double k = kvec.norm();

    // Analytic integration
    vec3cd rad = vec3cd::Zero();
    int triIdx = 0;

    for (const auto& tri : tris) {
    // const auto& tri = tris[triIdx];

        const auto& Xs = tri->Xs;
        const auto& Dps = tri->Dps;

        const auto& dX = X - Xs[0];

        const double alpha = -kvec.dot(Dps[0]);
        const double beta = kvec.dot(Dps[2]);

        const cmplx expI_alpha = exp(iu*alpha);
        const cmplx expI_beta = exp(iu*beta);

        const cmplx I1 = (1.0 - expI_alpha) / alpha - (1.0 - expI_beta) / beta;
        const cmplx Ipart = 1.0/(alpha-beta) * I1;
        const cmplx I2 = -iu * (-Ipart - (1.0 - (1.0 - iu*alpha) * expI_alpha) / (alpha*alpha) );
        const cmplx I3 = -iu * (Ipart + (1.0 - (1.0 - iu*beta) * expI_beta) / (beta*beta));

        rad += exp(iu*kvec.dot(dX)) / (alpha-beta) * (I1 * (Xs[0] - Xpm[triIdx]) + I2*Dps[0] - I3*Dps[2]);

        if (triIdx++) rad *= -1.0;
    }
    //

    // Numeric integration
    vec3cd radNum = vec3cd::Zero();
    triIdx = 0;

    for (const auto& tri : tris) {

        auto [nodes, weight] = tri->getQuads();
        for (const auto& node : nodes)
            radNum += weight * exp(iu*kvec.dot(X-node))
            * (node - Xpm[triIdx]);

        if (triIdx++) radNum *= -1.0;
    }
    //

    std::cout << std::setprecision(9) << rad << '\n' << radNum << '\n';

    return leng * rad;
}

/* getRadAtPoint(X)
 * Return the full radiated field with given wavenum
 * due to this RWG at X
 */
vec3cd RWG::getRadAtPoint(const vec3d& X) const {

    vec3cd rad = vec3cd::Zero();

    const double k = Einc->wavenum;
   
    /*auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus) {

        const auto& dyadic = Math::dyadicG(X - quadNode, wavenum);

        rad += weightPlus * dyadic * (vPlus - quadNode);
    }

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus) {
    
        const auto& dyadic = Math::dyadicG(X - quadNode, wavenum);
        
        rad += weightMinus * dyadic * (quadNode - vMinus);
    }*/

    return leng * rad;
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx RWG::getIntegratedRad(const std::shared_ptr<Source> src) const {

    cmplx integratedRad = 0.0;

    const double k = Einc->wavenum;

    /*auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        integratedRad += weightPlus * src->getRadAtPoint(quadNode).dot(vPlus - quadNode);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        integratedRad += weightMinus * src->getRadAtPoint(quadNode).dot(quadNode - vMinus);*/

    return leng * integratedRad;
    
}