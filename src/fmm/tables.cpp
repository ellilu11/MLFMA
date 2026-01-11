#include "tables.h"

realVec FMM::Tables::dists;
std::vector<vec3d> FMM::Tables::rhats;
std::array<vec3d,316> FMM::Tables::dXs;

void FMM::Tables::buildAngularTables() {
    const auto& angles = Node::angles[level];
    const auto [nth, nph] = angles.getNumAngles();
    const size_t nDirs = nth*nph;

    khat.resize(nDirs);
    toThPh.resize(nDirs);
    ImRR.resize(nDirs);

    size_t iDir = 0;
    for (int ith = 0; ith < nth; ++ith) {
        const double theta = angles.thetas[ith];

        for (int iph = 0; iph < nph; ++iph) {
            const double phi = angles.phis[iph];

            khat[iDir] = Math::fromSph(vec3d(1.0, theta, phi));
            toThPh[iDir] = Math::toThPh(theta, phi);
            ImRR[iDir] = Math::ImRR(khat[iDir]);

            ++iDir;
        }
    }
}

std::vector<interpPair> FMM::Tables::getInterpTheta(int srcLvl, int tgtLvl) 
{
    const int order = Node::config.interpOrder;

    const realVec& srcThetas = Node::angles[srcLvl].thetas, tgtThetas = Node::angles[tgtLvl].thetas;
    const int mth = srcThetas.size(), nth = tgtThetas.size();

    std::vector<interpPair> interpPairs;

    for (size_t jth = 0; jth < nth; ++jth) {
        const double tgtTheta = tgtThetas[jth];

        const int nearIdx = Interp::getNearGLNodeIdx(tgtTheta, mth, 0.0, PI);

        // Assemble source thetas interpolating target theta
        realVec interpThetas(2*order);
        for (int ith = nearIdx+1-order, k = 0; ith <= nearIdx+order; ++ith, ++k) {

            // Flip ith if not in [0, mth-1]
            int ith_flipped = Math::flipIdxToRange(ith, mth);

            double srcTheta = srcThetas[ith_flipped];

            // Extend source thetas to outside [0, pi] as needed
            if (ith < 0) srcTheta *= -1.0;
            else if (ith >= mth) srcTheta = 2.0*PI - srcTheta;

            interpThetas[k] = srcTheta;
        }

        vecXd coeffs(2*order);
        for (int k = 0; k < 2*order; ++k)
            coeffs[k] = 
                Interp::evalLagrangeBasis(tgtTheta, interpThetas, k);

        interpPairs.emplace_back(coeffs, nearIdx);

    }

    return interpPairs;
}

std::vector<interpPair> FMM::Tables::getInterpPhi(int srcLvl, int tgtLvl)
{
    const int order = Node::config.interpOrder;

    const realVec& srcPhis = Node::angles[srcLvl].phis, tgtPhis = Node::angles[tgtLvl].phis;
    const int mph = srcPhis.size(), nph = tgtPhis.size();

    std::vector<interpPair> interpPairs;

    for (size_t jph = 0; jph < nph; ++jph) {
        const double tgtPhi = tgtPhis[jph];

        const int nearIdx = std::floor(mph * tgtPhi / (2.0*PI));

        // Assemble source phis interpolating target phi
        realVec interpPhis(2*order);
        for (int iph = nearIdx+1-order, k = 0; iph <= nearIdx+order; ++iph, ++k)
            interpPhis[k] = 2.0*PI*iph/static_cast<double>(mph);

        vecXd coeffs(2*order);
        for (int k = 0; k < 2*order; ++k)
            coeffs[k] =
                Interp::evalLagrangeBasis(tgtPhi, interpPhis, k);

        interpPairs.emplace_back(coeffs, nearIdx);
    }

    return interpPairs;
}

void FMM::Tables::buildInterpTables() {
    // M2M interpolation tables
    interpTheta = getInterpTheta(level+1, level);
    interpPhi = getInterpPhi(level+1, level);

    // L2L interpolation tables
    //invInterpTheta = getInterpTheta(level, level+1);
    //invInterpPhi = getInterpPhi(level, level+1);
}

Map<vecXcd> FMM::Tables::getAlpha() {
    using namespace Math;

    const double wavenum = Node::wavenum;

    const int L = Node::angles[level].L;
    const int nth = Node::angles[level].getNumAngles().first;
    const int nps = std::floor(Node::config.overInterp*(nth-1));
    const double nodeLeng = Node::config.rootLeng / pow(2.0, level);

    Map<vecXcd> alpha;
        
    for (const auto& dist : dists) {
        const double kr = wavenum * dist * nodeLeng;

        vecXcd transl_dist(nps);

        for (int ips = 0; ips < nps; ++ips) {

            const double xi = cos(PI*ips/static_cast<double>(nps-1));
            cmplx coeff = 0.0;

            for (int l = 0; l <= L; ++l)
                coeff += 
                    powI(l) * (2.0*l+1.0)
                    * sphericalHankel1(kr, l) 
                    * legendreP(xi, l).first;

            transl_dist[ips] = iu * wavenum / (4.0*PI) * coeff;
        }

        alpha.emplace(dist, transl_dist);
    }

    assert(alpha.size() == dists.size());

    return alpha;
};

HashMap<interpPair> FMM::Tables::getInterpPsi() {
    const int order = Node::config.interpOrder;

    // Find all unique psi = acos(khat.dot(rhat))
    const auto [nth, nph] = Node::angles[level].getNumAngles();
    const int nDir = nth*nph;

    realVec psis(nDir*rhats.size());

    size_t m = 0;
    for (size_t iDir = 0; iDir < nDir; ++iDir) {
        const auto& khat = this->khat[iDir];

        for (const auto& rhat : rhats)
            psis[m++] = acos(khat.dot(rhat));
    }

    std::sort(psis.begin(), psis.end()); 
    psis.erase(std::unique(psis.begin(), psis.end()), psis.end());

    // Compute Lagrange coefficients for each possible psi
    const int nps = std::floor(Node::config.overInterp*(nth-1));

    HashMap<interpPair> interpPairs;

    size_t idx = 0;
    for (auto psi : psis) {
        // Find idx of psi node nearest this psi
        const int nearIdx = std::floor((nps-1) * psi / PI);

        // Assemble psis interpolating this psi
        realVec psis;
        for (int ips = nearIdx+1-order; ips <= nearIdx+order; ++ips)
            psis.push_back(PI*ips/static_cast<double>(nps-1));

        // CONSIDER: Use barycentric coordinates
        vecXd coeffs(2*order);
        for (size_t k = 0; k < 2*order; ++k)
            coeffs[k] = Interp::evalLagrangeBasis(psi, psis, k); 

        interpPairs.emplace(psi, std::make_pair(coeffs, nearIdx));

        // if (!level) std::cout << std::setprecision(15) << psi << '\n';
    }

    // assert(interpPairs.size() == psis.size());
    // std::cout << interpPairs.size() << ' ' << psis.size() << '\n';

    return interpPairs;
}

void FMM::Tables::buildTranslationTable() {
    const int order = Node::config.interpOrder;

    const auto& alphas = getAlpha();
    const auto& interpPsis = getInterpPsi(); 

    const auto [nth, nph] = Node::angles[level].getNumAngles();
    const int nDir = nth*nph;

    const int nps = std::floor(Node::config.overInterp*(nth-1));

    for (const auto& dX : dXs) {
        const double r = dX.norm();
        const auto& rhat = dX / r;

        const auto& alpha_dX = alphas.at(r);

        vecXcd transl_dX(nth*nph);

        for (int iDir = 0; iDir < nDir; ++iDir) {
            const auto& khat = this->khat[iDir];

            const double psi = acos(khat.dot(rhat));

            cmplx translCoeff = 0.0;

            const auto [interpPsi, nearIdx] = interpPsis.at(psi);

            for (int ips = nearIdx+1-order, k = 0; k < 2*order; ++ips, ++k) {

                const int ips_flipped = Math::flipIdxToRange(ips, nps);

                translCoeff += alpha_dX[ips_flipped] * interpPsi[k];
            }

            transl_dX[iDir] = translCoeff;
        }

        transl.emplace(dX, transl_dX);
    }

    assert(transl.size() == dXs.size());

}