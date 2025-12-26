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
   
    /* Analytic integration
    vec3cd rad = vec3cd::Zero();
    int triIdx = 0;

    for (const auto& tri : tris) {

        const auto& Xs = tri->Xs, Ds = tri->Ds;

        const auto& dX = X - Xs[0];

        const double 
            alpha = -kvec.dot(Ds[0]), 
            beta = kvec.dot(Ds[2]),
            gamma = alpha-beta;

        const double alphasq = alpha*alpha, betasq = beta*beta;

        // std::cout << alpha << ' ' << beta << ' ' << gamma << '\n';

        const cmplx expI_alpha = exp(iu*alpha), expI_beta = exp(iu*beta);

        const cmplx 
            f0_alpha = (alpha ? (1.0 - expI_alpha) / alpha : -iu), 
            f0_beta = (beta ? (1.0 - expI_beta) / beta : -iu);

        const cmplx 
            f1_alpha = (alpha ? (1.0 - (1.0 - iu*alpha) * expI_alpha) / alphasq : -0.5),
            f1_beta = (beta ? (1.0 - (1.0 - iu*beta) * expI_beta) / betasq : -0.5);

        const cmplx
            f2 = (expI_alpha*(alphasq + 2.0*iu*alpha - 2.0) + 2.0) / (2.0*pow(alpha,3)); // TODO: Debug

        const cmplx
            I1 = (gamma ? (f0_alpha - f0_beta) / gamma : -f1_alpha),
            I2 = -iu * (gamma ? (-I1 - f1_alpha) / gamma : -f2 ),
            I3 = -iu * (gamma ? (I1 + f1_beta) / gamma : f2 );

        rad += exp(iu*kvec.dot(dX))
                * (I1 * (Xs[0] - Xpm[triIdx]) + I2*Ds[0] - I3*Ds[2])
                * Math::sign(triIdx);

        ++triIdx;

        //if (!gamma) {
        //    // std::cout << Ds[0] << ' ' << Ds[2] << ' ' << I2*Ds[0] - I3*Ds[2] << '\n';
        //    std::cout << "Anl: " << std::setprecision(9) << rad << '\n';
        //}
    }
    */

    // Numeric integration
    vec3cd radNum = vec3cd::Zero();
    int triIdx = 0;

    for (const auto& tri : tris) {

        auto [nodes, weight] = tri->getQuads();
        for (const auto& node : nodes)
            radNum += weight * exp(iu*kvec.dot(X-node))
                        * (node - Xpm[triIdx])
                        * Math::sign(triIdx);
        ++triIdx;

        // if (!gamma) std::cout << "Num: " << std::setprecision(9) << radNum << '\n';
        
    }
    //

    // std::cout << std::setprecision(9) << rad << '\n' << radNum << "\n\n";
    // std::cout << std::setprecision(9) << (rad - radNum).norm() / radNum.norm() << '\n';

    return leng * radNum;
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx RWG::getIntegratedRad(const std::shared_ptr<Source> src) const {

    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    const double k = Einc->wavenum;

    cmplx intRad = 0.0;

    int obsTriIdx = 0;
    for (const auto& obsTri : tris) {

        auto [obsNodes, obsWeight] = obsTri->getQuads(); // TODO: Optimize
        const auto& obsXpm = Xpm[obsTriIdx];

        int srcTriIdx = 0;
        for (const auto& srcTri : srcRWG->tris) {

            if (obsTri == srcTri) continue; // TODO: Handle coincident tris

            // if (obsTri->isAdjacent(srcTri)) continue; // TODO: Handle adjacent tris

            auto [srcNodes, srcWeight] = srcTri->getQuads();
            const auto& srcXpm = srcRWG->Xpm[srcTriIdx];

            for (const auto& obs : obsNodes) {
                for (const auto& src : srcNodes) {

                    const vec3cd& rad = srcWeight 
                        * Math::dyadicG(obs-src, k) * (src - srcXpm) 
                        * Math::sign(srcTriIdx);

                    // intRad += obsWeight * (obs - obsXpm).dot(rad.conjugate()) // Hermitian dot!
                    intRad += obsWeight * conj(rad.dot(obs-obsXpm)) // Hermitian dot!
                        * Math::sign(obsTriIdx);

                }
            }

            ++srcTriIdx;
        }

        ++obsTriIdx;
    }

    return leng * srcRWG->leng * intRad;
}