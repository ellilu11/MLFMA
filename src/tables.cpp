#include "tables.h"

void Tables::buildAngularTables(
    const int maxLevel,
    const std::vector<realVec>& thetas,
    const std::vector<realVec>& phis,
    const double wavenum) {
    for (int level = 0; level <= maxLevel; ++level) {

        const int nth = thetas[level].size();
        const int nph = phis[level].size();

        std::vector<mat3d> ImKK_lvl;
        std::vector<vec3d> kvec_lvl;
        std::vector<Eigen::Matrix<double,2,3>> matToThPh_lvl;
        std::vector<Eigen::Matrix<double,3,2>> matFromThPh_lvl;
        
        std::vector<mat3d> matToSph_lvl;
        std::vector<mat3d> matFromSph_lvl;

        for (int ith = 0; ith < nth; ++ith) {
            const double th = thetas[level][ith];
            for (int iph = 0; iph < nph; ++iph) {
                const double ph = phis[level][iph];

                ImKK_lvl.push_back(Math::IminusRR(th, ph));
                kvec_lvl.push_back(Math::vecSph(wavenum, th, ph));
                matToThPh_lvl.push_back(Math::matToThPh(th, ph));
                matFromThPh_lvl.push_back(Math::matFromThPh(th, ph));

                matToSph_lvl.push_back(Math::matToSph(th, ph));
                matFromSph_lvl.push_back(Math::matFromSph(th, ph));
            }
        }

        ImKK.push_back(ImKK_lvl);
        kvec.push_back(kvec_lvl);
        matToThPh.push_back(matToThPh_lvl);
        matFromThPh.push_back(matFromThPh_lvl);

        matToSph.push_back(matToSph_lvl);
        matFromSph.push_back(matFromSph_lvl);
    }
}

void Tables::buildInterpThetaTable(
    const int maxLevel, const int order, const std::vector<realVec>& thetas) {

    // Do not construct interp table for leaf level
    for (size_t level = 0; level < maxLevel; ++level) { 
        std::vector<realVec> interpTheta_lvl;
        std::vector<int> t_lvl;

        const int nth = thetas[level].size();
        const int mth = thetas[level+1].size();

        for (size_t jth = 0; jth < nth; ++jth) {
            realVec interpTheta_lvl_th;
            const double theta = thetas[level][jth];

            // Find idx of child theta nearest parent theta
            const int t = Interp::getNearGLNodeIdx(theta, mth, 0.0, PI);

            // Assemble child thetas interpolating parent theta
            realVec branchThetas;
            for (int ith = t+1-order; ith <= t+order; ++ith) {

                // Flip ith if not in [0, mth-1]
                int ith_flipped = Math::flipIdxToRange(ith, mth);

                auto branchTheta = thetas[level+1][ith_flipped];

                // Extend interpolating thetas to outside [0, pi] as needed
                if (ith < 0)
                    branchTheta *= -1.0;
                else if (ith >= mth)
                    branchTheta = 2.0*PI - branchTheta;

                branchThetas.push_back(branchTheta);

                // std::cout << jth << ' ' << ith << ' ' << branchTheta << '\n';
            }

            // std::cout << '\n';

            for (int k = 0; k < 2*order; ++k) {
                interpTheta_lvl_th.push_back(
                    Interp::evalLagrangeBasis(theta, branchThetas, k));

                // std::cout << ith << ' ' << k << ' ' << interpTheta_lvl_th[k] << '\n';
            }

            // std::cout << '\n';

            interpTheta_lvl.push_back(interpTheta_lvl_th);
            t_lvl.push_back(t);
        }

        interpTheta.push_back(interpTheta_lvl);

        ts.push_back(t_lvl);
    }
}

void Tables::buildInterpPhiTable(
    const int maxLevel, const int order, const std::vector<realVec>& phis) {

    // Do not construct interp table for leaf level
    for (size_t level = 0; level < maxLevel; ++level) {
        std::vector<realVec> interpPhi_lvl;
        std::vector<int> s_lvl;

        const int nph = phis[level].size();
        const int mph = phis[level+1].size();

        for (size_t jph = 0; jph < nph; ++jph) {
            realVec interpPhi_lvl_ph;
            const double phi = phis[level][jph];

            // Find idx of child phi nearest parent phi
            const int s = std::floor(mph * phi / (2.0 * PI));

            // Assemble child phis interpolating parent phi
            realVec branchPhis;
            for (int iph = s+1-order; iph <= s+order; ++iph) {

                // Wrap iph if not in [0, mph-1]
                // size_t iph_wrapped = Math::wrapIdxToRange(iph, mph); 
                // auto branchPhi = phis[level+1][iph_wrapped];

                auto branchPhi = 2.0*PI*iph/static_cast<double>(mph);

                branchPhis.push_back(branchPhi);
            }

            for (size_t k = 0; k < 2*order; ++k)
                interpPhi_lvl_ph.push_back(
                    Interp::evalLagrangeBasis(phi, branchPhis, k));
            
            interpPhi_lvl.push_back(interpPhi_lvl_ph);
            s_lvl.push_back(s);
        }

        interpPhi.push_back(interpPhi_lvl);
        ss.push_back(s_lvl);
    }
}

void Tables::buildTranslationTable(
    const int maxLevel, const int order, const std::vector<int>& Ls, 
    const double rootLeng, const double wavenum) {

    using namespace Math;

    const auto& dists = getINodeDistances();

    for (size_t level = 0; level <= maxLevel; ++level) {

        const auto L = Ls[level];
        const int nps = std::floor(q*L);
        const double nodeLeng = rootLeng / pow(2.0, level);

        std::map<double,cmplxVec,Comp> transl_lvl;
        
        for (const auto& dist : dists) {

            cmplxVec transl_dist;

            for (int ips = 0; ips < nps; ++ips) {
                const double psi = PI * ips / static_cast<double>(nps-1);
                cmplx coeff = 0.0;

                for (int l = 0; l <= L; ++l) {

                    coeff += powI(l) * (2.0*l+1.0)
                        * sphericalHankel1(wavenum*dist*nodeLeng, l)
                        * legendreP(cos(psi), l).first
                        ;

                    /*if (level == 0 && iDist == 0 && ips == 0)
                        std::cout << l << ' ' << ' ' << 
                        term << '\n';
                     */
                }

                transl_dist.push_back(coeff);
            }

            transl_lvl.emplace(dist, transl_dist);
        }

        transl.push_back(transl_lvl);
    }

};

void Tables::buildInterpPsiTable(
    const std::vector<realVec>& thetas,
    const std::vector<realVec>& phis,
    int maxLevel,
    int order, 
    double wavenum) {

    const auto& dirs = Math::getINodeDirections();

    for (size_t level = 0; level <= maxLevel; ++level) {
        const int nth = thetas[level].size();
        const int nph = phis[level].size();

        const auto L = nth - 1;
        const int nps = std::floor(q*L);

        std::map<double, realVec, Comp> interpPsi_lvl;

        size_t idx = 0;
        for (size_t ith = 0; ith < nth; ++ith) {
            for (size_t iph = 0; iph < nph; ++iph) {

                const auto& khat = kvec[level][idx] / wavenum;

                // Loop over all possible rhat
                for (const auto& rhat : dirs) {

                    const double psi = acos(khat.dot(rhat));

                    // Find idx of psi node nearest this psi
                    const int s = std::floor(nps * psi / PI);

                    // Assemble psis interpolating this psi
                    realVec psis;
                    for (int ips = s+1-order; ips <= s+order; ++ips)
                        psis.push_back(PI*ips/static_cast<double>(nps));

                    realVec interpPsi_ps;

                    for (size_t k = 0; k <= 2*order-1; ++k)
                        interpPsi_ps.push_back(
                            Interp::evalLagrangeBasis(psi, psis, k));

                    interpPsi_lvl.emplace(psi, interpPsi_ps);

                }
            }
        }

        interpPsi.push_back(interpPsi_lvl);
    }

}