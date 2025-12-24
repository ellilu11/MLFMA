#include "tables.h"

void Tables::buildAngularTables() {
   
    for (int level = 0; level <= Node::maxLevel; ++level) {

        const auto [nth, nph] = Node::getNumAngles(level);

        std::vector<vec3d> khat_lvl(nth*nph);
        std::vector<mat23d> toThPh_lvl(nth*nph);

        size_t idx = 0;
        for (int ith = 0; ith < nth; ++ith) {
            const double th = Node::thetas[level][ith];

            for (int iph = 0; iph < nph; ++iph) {
                const double ph = Node::phis[level][iph];

                khat_lvl[idx] = Math::fromSph(vec3d(1.0, th, ph));
                toThPh_lvl[idx++] = Math::toThPh(th, ph);
            }
        }

        khat.push_back(khat_lvl);
        toThPh.push_back(toThPh_lvl);
    }
}

std::vector<interpPair> Tables::getInterpThetaAtLvl(int srcLvl, int tgtLvl) {

    assert(abs(srcLvl - tgtLvl) == 1);

    const int mth = Node::getNumAngles(srcLvl).first;
    const int nth = Node::getNumAngles(tgtLvl).first;

    const auto& srcThetas = Node::thetas[srcLvl];
    const auto& tgtThetas = Node::thetas[tgtLvl];

    std::vector<interpPair> interpPairs;

    for (size_t jth = 0; jth < nth; ++jth) {
        const double tgtTheta = tgtThetas[jth];

        const int nearIdx = Interp::getNearGLNodeIdx(tgtTheta, mth, 0.0, PI);

        // Assemble source thetas interpolating target theta
        realVec interpThetas;
        for (int ith = nearIdx+1-order; ith <= nearIdx+order; ++ith) {

            // Flip ith if not in [0, mth-1]
            int ith_flipped = Math::flipIdxToRange(ith, mth);

            double srcTheta = srcThetas[ith_flipped];

            // Extend source thetas to outside [0, pi] as needed
            if (ith < 0)
                srcTheta *= -1.0;
            else if (ith >= mth)
                srcTheta = 2.0*PI - srcTheta;

            interpThetas.push_back(srcTheta);
        }

        vecXd coeffs(2*order);
        for (int k = 0; k < 2*order; ++k)
            coeffs[k] = 
                Interp::evalLagrangeBasis(tgtTheta, interpThetas, k);

        interpPairs.emplace_back(coeffs, nearIdx);

    }

    return interpPairs;
}

std::vector<interpPair> Tables::getInterpPhiAtLvl(int srcLvl, int tgtLvl) {

    assert(abs(srcLvl - tgtLvl) == 1);

    const int mph = Node::getNumAngles(srcLvl).second;
    const int nph = Node::getNumAngles(tgtLvl).second;

    const auto& srcPhis = Node::phis[srcLvl];
    const auto& tgtPhis = Node::phis[tgtLvl];

    std::vector<interpPair> interpPairs;

    for (size_t jph = 0; jph < nph; ++jph) {
        const double tgtPhi = tgtPhis[jph];

        const int nearIdx = std::floor(mph * tgtPhi / (2.0*PI));

        // Assemble source phis interpolating target phi
        realVec interpPhis;
        for (int iph = nearIdx+1-order; iph <= nearIdx+order; ++iph)
            interpPhis.push_back(2.0*PI*iph/static_cast<double>(mph));

        vecXd coeffs(2*order);
        for (int k = 0; k < 2*order; ++k)
            coeffs[k] = 
                Interp::evalLagrangeBasis(tgtPhi, interpPhis, k);

        interpPairs.emplace_back(coeffs, nearIdx);
    }

    return interpPairs;
}

void Tables::buildInterpTables() {

    for (int lvl = 0; lvl < Node::maxLevel; ++lvl) {

        // M2M interpolation tables
        interpTheta.push_back(getInterpThetaAtLvl(lvl+1, lvl));
        interpPhi.push_back(getInterpPhiAtLvl(lvl+1, lvl));

        // L2L interpolation tables
        invInterpTheta.push_back(getInterpThetaAtLvl(lvl, lvl+1));
        invInterpPhi.push_back(getInterpPhiAtLvl(lvl, lvl+1));

    }
}

void Tables::buildTranslationTable() {

    using namespace Math;

    const int rootLeng = Node::config.rootLeng;
    const double wavenum = Node::wavenum;

    const auto& dists = getINodeDistances();

    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const int L = Node::Ls[level];

        const int nth = Node::getNumAngles(level).first;
        const int nps = std::floor(Node::config.overInterp*(nth-1));

        const double nodeLeng = rootLeng / pow(2.0, level);

        Map<double,vecXcd> transl_lvl;
        
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

            transl_lvl.emplace(dist, transl_dist);
        }

        transl.push_back(transl_lvl);
    }

};

void Tables::buildInterpPsiTable() {

    // auto start = Clock::now();

    const auto& rhats = Math::getINodeDirections();

    std::vector<realVec> psis(Node::maxLevel+1);

    // Find all unique psi = acos(khat.dot(rhat)) at each level
    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const auto [nth, nph] = Node::getNumAngles(level);

        realVec psis_lvl(nth*nph*rhats.size());

        size_t m = 0;
        for (size_t l = 0; l < nth*nph; ++l) {

            const auto& khat_ = khat[level][l];

            for (const auto& rhat : rhats)
                psis_lvl[m++] = acos(khat_.dot(rhat));
        }

        std::sort(psis_lvl.begin(), psis_lvl.end()); 
        psis_lvl.erase(
            std::unique(psis_lvl.begin(), psis_lvl.end()), psis_lvl.end());

        psis[level] = psis_lvl;
    }

    //auto end = Clock::now();
    //Time duration_ms = end - start;
    //std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // start = Clock::now();

    // Compute Lagrange coefficients for each possible psi at each level
    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const int nth = Node::getNumAngles(level).first;
        const int nps = std::floor(Node::config.overInterp*(nth-1));

        HashMap<double,interpPair> interpPairs;

        size_t idx = 0;
        for (auto psi : psis[level]) {

            // Find idx of psi node nearest this psi
            const int nearIdx = std::floor((nps-1) * psi / PI);

            // Assemble psis interpolating this psi
            realVec psis;
            for (int ips = nearIdx+1-order; ips <= nearIdx+order; ++ips)
                psis.push_back(PI*ips/static_cast<double>(nps-1));

            vecXd coeffs(2*order);
            for (size_t k = 0; k < 2*order; ++k)
                coeffs[k] = 
                     Interp::evalLagrangeBasis(psi, psis, k); // TODO: Use barycentric coordinates

            interpPairs.emplace(psi, std::make_pair(coeffs, nearIdx));

        }

        interpPsi.push_back(interpPairs);
    }

    //auto end = Clock::now();
    //auto duration_ms = end - start;
    //std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

}
