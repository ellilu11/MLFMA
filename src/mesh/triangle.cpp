#include "triangle.h"

std::vector<Mesh::quadPair> Mesh::Triangle::quadCoeffs;
int Mesh::Triangle::numQuads;

/* Triangle(iVerts)
 * Construct triangle from global indices of vertices
 * iVerts : global indices of vertices
 */
Mesh::Triangle::Triangle(const vec3i& iVerts, int iTri)
    : iVerts(iVerts), iTri(iTri)
{
    auto Xs = getVerts();
    std::array<vec3d, 3> Ds; // edge vectors
    for (int i = 0; i < 3; ++i) Ds[i] = Xs[(i+1)%3] - Xs[i];

    center = (Xs[0] + Xs[1] + Xs[2]) / 3.0;
    nhat = (Ds[0].cross(Ds[1])).normalized();
    area = (Ds[0].cross(-Ds[2])).norm() / 2.0;
    assert(area > 0.0);

    buildTriQuads();
    reverseOrient();
    // buildSelfIntegratedInvR();
}

/* buildTriQuads()
 * Build quadrature nodes and weights for this triangle
 * Use barycentric coordinates to transform from reference triangle to this triangle
 */
void Mesh::Triangle::buildTriQuads() {
    auto [X0, X1, X2] = getVerts();

    auto baryToPos = [&](const vec3d& ws) {
        return ws[0]*X0 + ws[1]*X1 + ws[2]*X2;
        };

    quads.reserve(quadCoeffs.size());
    for (const auto& [coeffs, weight] : quadCoeffs)
        quads.emplace_back(baryToPos(coeffs), weight);
}

/* reverseOrient()
 * If nhat is pointing inward, reverse it
 * Assume a closed, star-shaped mesh centered at and enclosing rootCenter
 */
void Mesh::Triangle::reverseOrient() {
    if ((center-rootCenter).dot(nhat) < 0.0) {
        // std::cout << "   Reversing normals of triangle " << iTri << '\n';
        nhat *= -1.0;
        std::swap(this->iVerts[0], this->iVerts[2]); // Swap verts per right-hand rule
    }
}

/* buildSelfIntegratedInvR()
 * Build self-integrated 1/R and its related integrals for this triangle
 */
void Mesh::Triangle::buildSelfIntegratedInvR() {
    auto [V0, V1, V2] = getVerts();
    double a00 = V0.dot(V0), a01 = V0.dot(V1), a02 = V0.dot(V2);
    double a11 = V1.dot(V1), a12 = V1.dot(V2), a22 = V2.dot(V2);
    double l0 = (V1-V2).norm(), l1 = (V2-V0).norm(), l2 = (V0-V1).norm();

    auto logN = [](double l0, double l1, double l2) {
        double lsum = l0+l1, ldiff = l2-l0;
        return std::log((lsum*lsum - l2*l2) / (l1*l1 - ldiff*ldiff));
        };
    double log0 = logN(l0, l1, l2), log1 = logN(l1, l2, l0), log2 = logN(l2, l0, l1);

    // Integral of \lambda_i * \lambda_i' / r
    auto f0 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log0 / (20.0*l0)
            + log1 * (l0sq + 5.0*l1sq - l2sq) / (120.0*l1sq*l1)
            + log2 * (l0sq - l1sq + 5.0*l2sq) / (120.0*l2sq*l2)
            + (l2-l0) / (60.0*l1sq) + (l1-l0) / (60.0*l2sq);
        };

    // Integral of \lambda_i * \lambda_j' / r
    auto f1 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log2 / (40.0*l2)
            + log0 * (3.0*l0sq + l1sq - l2sq) / (80.0*l0sq*l0)
            + log1 * (l0sq + 3.0*l1sq - l2sq) / (80.0*l1sq*l1)
            + (l2-l1) / (40.0*l0sq) + (l2-l0) / (40.0*l1sq);
        };

    // Integral of \lambda_i / r or \lambda_i' / r
    auto f2 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log0 / (8.0*l0)
            + log1 * (l0sq + 5.0*l1sq - l2sq) / (48.0*l1sq*l1)
            + log2 * (l0sq - l1sq + 5.0*l2sq) / (48.0*l2sq*l2)
            + (l2-l0) / (24.0*l1sq) + (l1-l0) / (24.0*l2sq);
        };

    selfInts[0]
        = f0(l0, l1, l2, log0, log1, log2) * (a00 - 2.0*a01 + a11)
        + f0(l1, l2, l0, log1, log2, log0) * (a00 - 2.0*a02 + a22)
        + f1(l0, l1, l2, log0, log1, log2) * (a00 - a02 - a01 + a12)
        + f1(l0, l1, l2, log0, log1, log2) * (a00 - a01 - a02 + a12);
    selfInts[1] = f2(l0, l1, l2, log0, log1, log2);
    selfInts[2] = f2(l1, l2, l0, log1, log2, log0);
    selfInts[3] = (log0/l0 + log1/l1 + log2/l2) / 3.0;
}

/* getDoubleSelfIntegratedInvR(vobs, vsrc)
 * Get the double integral of 1/|vobs-r'| over this triangle
 */
double Mesh::Triangle::getDoubleSelfIntegratedInvR(const vec3d& vobs, const vec3d& vsrc) const
{
    double k2 = config.k * config.k;
    const auto [V0, V1, V2] = getVerts();
    double a00 = V0.dot(V0), a01 = V0.dot(V1), a02 = V0.dot(V2); // cache?
    const vec3d& vsum = vsrc + vobs;
    
    return (selfInts[0]
        + selfInts[1] * (-2.0*a00 + 2.0*a01 + (V0-V1).dot(vsum))
        + selfInts[2] * (-2.0*a00 + 2.0*a02 + (V0-V2).dot(vsum))
        + selfInts[3] * (a00 - V0.dot(vsum) + vsrc.dot(vobs) - 4.0/k2)) / (4.0*PI);
}

/* getIntegratedPlaneWave(kvec)
 * Get the integral of exp(i kvec . r') dr' over this triangle
 * Return a pair of scalar and vector parts of the integral, 
 * where the scalar part is the integral of exp(i kvec . r') 
 * and the vector part is the integral of (r' - center) * exp(i kvec . r')
 */
std::pair<cmplx, vec3cd>
Mesh::Triangle::getIntegratedPlaneWave(const vec3d& kvec) const {
    auto [X0, X1, X2] = getVerts();
    vec3d D0 = X1-X0, D2 = X0-X2;

    double alpha = kvec.dot(D0), beta = -kvec.dot(D2), gamma = alpha-beta;
    double alphasq = alpha*alpha, betasq = beta*beta;
    cmplx expI_alpha = exp(iu*alpha), expI_beta = exp(iu*beta);
    cmplx // TODO: Only compute if gamma != 0
        f0_alpha = (lzero(alpha) ? -iu : (1.0 - expI_alpha) / alpha),
        f0_beta = (lzero(beta) ? -iu : (1.0 - expI_beta) / beta);
    cmplx
        f1_alpha = (lzero(alpha) ? -0.5 : (1.0 - (1.0 - iu*alpha) * expI_alpha) / alphasq),
        f1_beta = (lzero(beta) ? -0.5 : (1.0 - (1.0 - iu*beta) * expI_beta) / betasq);

    cmplx scaRad; 
    vec3cd vecRad;
    if (lzero(gamma)) {
        cmplx
            f2 = (lzero(alpha) ? iu/6.0 :
                (expI_alpha*(alphasq + 2.0*iu*alpha - 2.0) + 2.0) / (2.0*alpha*alphasq));
        scaRad = -f1_alpha;
        vecRad = -iu*f2 * (D0 - D2);
        // radVec = -f1_alpha * (Xs[0] - Xnc[iTri]) - iu*f2 * (Ds[0] - Ds[2]);
    } else {
        cmplx
            I0 = (f0_alpha - f0_beta) / gamma,
            I1 = iu * (I0 + f1_alpha),
            I2 = -iu * (I0 + f1_beta);
        scaRad = I0;
        vecRad = (I1*D0 - I2*D2) / gamma;
        // radVec = I0 * (Xs[0] - Xnc[iTri]) + (I1*Ds[0] - I2*Ds[2]) / gamma;
    }

    return std::make_pair(scaRad, vecRad);
}

/* getIntegratedInvR(obs, doNumeric)
 * Get the integral of 1/|obs-r'| dr' over this triangle
 * Return a pair of scalar and vector parts of the integral, 
 * where the scalar part is the integral of 1/|obs-r'| 
 * and the vector part is the integral of (proj(r') - proj(obs)) / |obs-r'|
 * where proj() is the projection onto the plane of the triangle
 * If doNumeric is true, use numerical quadrature instead of analytical formula
 */
std::pair<double, vec3d>
Mesh::Triangle::getIntegratedInvR(const vec3d& obs, bool doNumeric) const
{
    double scaRad = 0.0;
    vec3d vecRad = vec3d::Zero();
    const vec3d& R = proj(obs);

    if (doNumeric) {
        for (const auto& [node, weight] : quads) {
            const double dr = (obs-node).norm();
            scaRad += weight / dr;
            vecRad -= weight * (R-proj(node)) / dr;
        }

        return std::make_pair(scaRad, vecRad);
    }

    auto [X0, X1, X2] = getVerts();
    double d = std::fabs(nhat.dot(obs-X0)), dsq = d*d;
    std::array<vec3d, 3> Ps =
        { proj(X0)-R, proj(X1)-R, proj(X2)-R };

    for (int i = 0; i < 3; ++i) {
        const vec3d& P0 = Ps[i];
        const vec3d& P1 = Ps[(i+1)%3];

        const vec3d& lhat = (P1-P0).normalized();
        const vec3d& uhat = lhat.cross(nhat);
        double l0 = P0.dot(lhat), l1 = P1.dot(lhat);

        const vec3d& P = P0 - l0*lhat;
        double p = P.norm(), p0 = P0.norm(), p1 = P1.norm();
        const vec3d& phat = (lzero(p) ? vec3d::Zero() : P.normalized());

        double rsq = p*p + dsq, r0 = std::sqrt(p0*p0 + dsq), r1 = std::sqrt(p1*p1 + dsq);
  
        double f2 = 
            (lzero(r0+l0) || lzero(r1+l1) ?  // use && ?
                std::log(l0/l1) :
                std::log((r1+l1)/(r0+l0)));

        scaRad += phat.dot(uhat) * (
            p*f2 - d * (std::atan2(p*l1, rsq+d*r1) - std::atan2(p*l0, rsq+d*r0)));
        vecRad += uhat * (rsq*f2 + l1*r1 - l0*r0);
    }
    scaRad /= 2.0*area;
    vecRad /= 4.0*area;

    return std::make_pair(scaRad, vecRad);
}

/* getIntegratedInvRcubed(obs, doNumeric)
 * Get the integral of 1/|obs-r'|^3 dr' over this triangle
 * Return a pair of scalar and vector parts of the integral, 
 * where the scalar part is the integral of 1/|obs-r'|^3
 * and the vector part is the integral of (proj(r') - proj(obs)) / |obs-r'|^3
 * where proj() is the projection onto the plane of the triangle
 * If doNumeric is true, use numerical quadrature instead of analytical formula
 */
std::pair<double, vec3d>
Mesh::Triangle::getIntegratedInvRcubed(const vec3d& obs, bool doNumeric) const
{
    double scaRad3 = 0.0;
    vec3d vecRad3 = vec3d::Zero();
    const vec3d& obsProj = proj(obs);

    if (doNumeric) {
        for (const auto& [node, weight] : quads) {
            double dr3 = pow((obs-node).norm(), 3);
            scaRad3 += weight / dr3;
            vecRad3 -= weight * (obsProj-proj(node)) / dr3;
        }
        return std::make_pair(scaRad3, vecRad3);
    }

    const auto& Xs = getVerts();
    double d = std::fabs(nhat.dot(obs-Xs[0])), dsq = d*d;

    std::array<vec3d, 3> Ps =
        { proj(Xs[0])-obsProj, proj(Xs[1])-obsProj, proj(Xs[2])-obsProj };

    for (int i = 0; i < 3; ++i) {
        const vec3d& P0 = Ps[i];
        const vec3d& P1 = Ps[(i+1)%3];

        const vec3d& lhat = (P1-P0).normalized();
        const vec3d& uhat = lhat.cross(nhat);
        double l0 = P0.dot(lhat), l1 = P1.dot(lhat);

        const vec3d& P = P0 - l0*lhat;
        double p = P.norm(), p0 = P0.norm(), p1 = P1.norm();

        const vec3d& phat = (lzero(p) ? vec3d::Zero() : P.normalized());

        double rsq = p*p + dsq, r0 = std::sqrt(p0*p0 + dsq), r1 = std::sqrt(p1*p1 + dsq);

        scaRad3 -= phat.dot(uhat) *
            (lzero(d) ?
                (lzero(p) ? 0.0 : l1/(p*r1) - l0/(p*r0)) :
                (atan2(d*l1, p*r1) - atan2(d*l0, p*r0) + atan2(l0, p) - atan2(l1, p)) / d);
        vecRad3 -= uhat *
            (lzero(r0+l0) || lzero(r1+l1) ?  // use && ?
                std::log(l0/l1) :
                std::log((r1+l1)/(r0+l0)));
    }
    scaRad3 /= 2.0*area;
    vecRad3 /= 2.0*area;

    return std::make_pair(scaRad3, vecRad3);
}

/* getSingularEFIE(srcTri, triPair, vobs, vsrc)
 * Get the singular part of the EFIE integral
 * Use the precomputed integratedInvR from triPair to avoid numerical quadrature
 */
double Mesh::Triangle::getSingularEFIE(
    const Triangle& srcTri, const vec3d& vobs, const vec3d& vsrc, size_t iPair) const
{
    double k2 = config.k * config.k;
    vec3d vsrcProj = srcTri.proj(vsrc);
    const intRads& intInvR =
        (iTri <= srcTri.iTri ? glTriPairs.intsInvR[iPair] : glTriPairs.intsInvR2[iPair]);

    double rad = 0.0;
    size_t iObs = 0;
    for (const auto& [obs, obsWeight] : quads) {
        vec3d obsProj = srcTri.proj(obs);
        auto [scaInvR, vecInvR] = intInvR[iObs];

        rad += ((obs-vobs).dot(vecInvR+(obsProj-vsrcProj)*scaInvR) - 4.0/k2*scaInvR)
                * obsWeight;
        ++iObs;
    }

    return rad / (4.0*PI); // apply factor of 1/(4pi)
}

/* getSingularMFIE(srcTri, triPair, vobs, vsrc)
 * Get the singular part of the MFIE integral
 * Use the precomputed integratedInvR and integratedInvRcubed from triPair to avoid numerical quadrature
 */
double Mesh::Triangle::getSingularMFIE(
    const Triangle& srcTri, const vec3d& vobs, const vec3d& vsrc, size_t iPair) const
{
    double k2 = config.k * config.k;
    vec3d vsrcProj = srcTri.proj(vsrc);

    const intRads& intInvR =
        (iTri <= srcTri.iTri ? glTriPairs.intsInvR[iPair] : glTriPairs.intsInvR2[iPair]);
    const intRads& intInvRcubed =
        (iTri <= srcTri.iTri ? glTriPairs.intsInvRcubed[iPair] : glTriPairs.intsInvRcubed2[iPair]);

    double rad = 0.0;
    size_t iObs = 0;
    for (const auto& [obs, obsWeight] : quads) {
        // (obs-src) x (src-vsrc) = [R - (src-vsrc)] x (src-vsrc) = R x (src-vsrc)
        vec3d R = obs - vsrc;
        vec3d obsProj = srcTri.proj(obs);

        auto [scaInvR, vecInvR] = intInvR[iObs];
        auto [scaInvR3, vecInvR3] = intInvRcubed[iObs];

        vec3d sumInvR = 0.5*k2 * (vecInvR + (obsProj-vsrcProj)*scaInvR);
        vec3d sumInvR3 = vecInvR3 + (obsProj-vsrcProj)*scaInvR3;

        rad += (obs-vobs).dot(nhat.cross(R.cross(sumInvR + sumInvR3)))
            * obsWeight;

        ++iObs;
    }

    return rad / (4.0*PI); // apply factor of 1/(4pi)
}