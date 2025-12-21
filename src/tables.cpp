#include "tables.h"

void Tables::buildAngularTables() {
   
    for (int level = 0; level <= Node::maxLevel; ++level) {

        const auto [nth, nph] = Node::getNumAngles(level);

        std::vector<vec3d> khat_lvl;
        std::vector<mat23d> toSphK_lvl;
        std::vector<mat23d> toSphKK_lvl;

        for (int ith = 0; ith < nth; ++ith) {
            const double th = Node::thetas[level][ith];

            for (int iph = 0; iph < nph; ++iph) {
                const double ph = Node::phis[level][iph];

                khat_lvl.push_back(Math::fromSph(vec3d(1.0, th, ph)));
                toSphK_lvl.push_back(Math::toSphR(th, ph));
                toSphKK_lvl.push_back(Math::toSphRR(th, ph));
            }
        }

        khat.push_back(khat_lvl);
        toSphK.push_back(toSphK_lvl);
        toSphKK.push_back(toSphK_lvl);
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
        realVec coeffs(2*order);
        const double tgtTheta = tgtThetas[jth];

        const int inear = Interp::getNearGLNodeIdx(tgtTheta, mth, 0.0, PI);

        // Assemble source thetas interpolating target theta
        // TODO: Use splicing with std::span
        realVec interpThetas;
        for (int ith = inear+1-order; ith <= inear+order; ++ith) {

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

        for (int k = 0; k < 2*order; ++k)
            coeffs[k] = 
                Interp::evalLagrangeBasis(tgtTheta, interpThetas, k);

        interpPairs.emplace_back(coeffs, inear);

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
        realVec coeffs(2*order);
        const double tgtPhi = tgtPhis[jph];

        const int inear = std::floor(mph * tgtPhi / (2.0*PI));

        // Assemble source phis interpolating target phi
        // TODO: Use splicing with std::span
        realVec interpPhis;
        for (int iph = inear+1-order; iph <= inear+order; ++iph)
            interpPhis.push_back(2.0*PI*iph/static_cast<double>(mph));

        for (int k = 0; k < 2*order; ++k)
            coeffs[k] = 
                Interp::evalLagrangeBasis(tgtPhi, interpPhis, k);

        interpPairs.emplace_back(coeffs, inear);
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

                const double psi = PI*ips/static_cast<double>(nps-1);
                const double xi = cos(psi);
                // const double xi = 2.0*ips/static_cast<double>(nps-1)-1.0;
                // const double xi = -2.0*ips/static_cast<double>(nps-1)+1.0;
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


void Tables::buildInterpPsiTable() { // CONSIDER: Interp over xi = cos(psi)

    //std::cout << "   Finding all psi\n";
    //auto start = Clock::now();

    const auto& rhats = Math::getINodeDirections();

    std::vector<realVec> psis;

    // Find all unique psi = acos(khat.dot(rhat)) at each level
    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const auto [nth, nph] = Node::getNumAngles(level);

        realVec psis_lvl;

        size_t idx = 0;
        for (size_t ith = 0; ith < nth; ++ith) {
            for (size_t iph = 0; iph < nph; ++iph) {

                const auto& khat_ = khat[level][idx++];

                // Loop over all possible rhat
                for (const auto& rhat : rhats)
                    psis_lvl.push_back(acos(khat_.dot(rhat)));
                    // psis_lvl.push_back(khat_.dot(rhat));

            }
        }

        std::sort(psis_lvl.begin(), psis_lvl.end()); 

        psis_lvl.erase(
            std::unique(psis_lvl.begin(), psis_lvl.end()), psis_lvl.end());

        // for (auto psi : psis_lvl) std::cout << psi << ' ';

        psis.push_back(psis_lvl);
    }

    //auto end = Clock::now();
    //Time duration_ms = end - start;
    //std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    //std::cout << "   Computing Lagrange coeffs for each psi\n";
    //start = Clock::now();

    // Compute Lagrange coefficients for each possible psi at each level
    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const int nth = Node::getNumAngles(level).first;
        const int nps = std::floor(Node::config.overInterp*(nth-1));

        HashMap<double,vecXcd> interpPsi_lvl;
        HashMap<double,int> s_lvl;

        size_t idx = 0;
        for (auto psi : psis[level]) {

            // Find idx of psi node nearest this psi
            const int s = std::floor((nps-1) * psi / PI);
            // const int s = std::floor((nps-1) * (psi + 1.0) / 2.0);
            // const int s = std::floor(-(nps-1) * (psi - 1.0) / 2.0);

            // Assemble psis interpolating this psi
            // TODO: Use splicing with std::span
            realVec psis;
            for (int ips = s+1-order; ips <= s+order; ++ips)
                psis.push_back(PI*ips/static_cast<double>(nps-1));
                // psis.push_back(2.0*ips/static_cast<double>(nps-1)-1.0);
                // psis.push_back(-2.0*ips/static_cast<double>(nps-1)+1.0);

            vecXd interpPsi_ps(2*order);

            for (size_t k = 0; k < 2*order; ++k)
                interpPsi_ps[k] = 
                     Interp::evalLagrangeBasis(psi, psis, k);

            interpPsi_lvl.emplace(psi, interpPsi_ps);
            s_lvl.emplace(psi, s);

        }

        interpPsi.push_back(interpPsi_lvl);
        idxPsi.push_back(s_lvl);
    }

    //end = Clock::now();
    //duration_ms = end - start;
    //std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

}
