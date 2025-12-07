#include "tables.h"

void Tables::buildAngularTables(
    const int maxLevel,
    const std::vector<realVec>& thetas,
    const std::vector<realVec>& phis,
    const double wavenum) {
    for (int level = 0; level <= maxLevel; ++level) {

        const int nth = thetas[level].size();
        const int nph = 2*nth;

        std::vector<mat3d> ImKK_lvl;
        std::vector<vec3d> kvec_lvl;
        std::vector<mat23d> matToThPh_lvl;
        // std::vector<mat3d> matToThPh_lvl;

        for (int ith = 0; ith < nth; ++ith) {
            const double th = thetas[level][ith];
            for (int iph = 0; iph < nph; ++iph) {
                const double ph = phis[level][iph];

                ImKK_lvl.push_back(Math::IminusRR(th, ph));
                kvec_lvl.push_back(Math::vecSph(wavenum, th, ph));
                matToThPh_lvl.push_back(Math::matToThPh(th, ph));
                // matToThPh_lvl.push_back(Math::matToSph(th, ph));
            }
        }

        ImKK.push_back(ImKK_lvl);
        kvec.push_back(kvec_lvl);
        matToThPh.push_back(matToThPh_lvl);
    }
}

void Tables::buildInterpThetaTable(
    const int maxLevel, const int order, const std::vector<realVec>& thetas) {

    // Do not construct interp table for leaf level
    for (size_t level = 0; level < maxLevel; ++level) { 
        std::vector<realVec> interpTheta_lvl;
        std::vector<size_t> t_lvl;

        const int nth = thetas[level].size();
        const int mth = thetas[level+1].size();

        for (size_t ith = 0; ith < nth; ++ith) {
            realVec interpTheta_lvl_th;
            const double theta = thetas[level][ith];

            // Find idx of child theta nearest parent theta
            const int t = Interp::getNearGLNodeIdx(theta, mth, 0.0, PI);

            // if (level) std::cout << ith << ' ' << theta << ' ' << t << '\n';

            // Assemble child thetas interpolating parent theta
            realVec branchThetas;
            for (int jth = t+1-order; jth <= t+order; ++jth) {

                // Flip jth if not in [0, mth-1]
                size_t jth_flipped = Math::flipIdxToRange(jth, mth);
                auto branchTheta = thetas[level+1][jth_flipped];

                // Extend interpolating thetas to outside [0, \pi] as needed
                if (jth < 0)
                    branchTheta *= -1.0;
                else if (jth >= mth)
                    branchTheta = 2.0*PI - branchTheta;

                branchThetas.push_back(branchTheta);

                // if (level) std::cout << ith << ' ' << jth << ' ' << branchTheta << '\n';
            }
            // if (level) std::cout << '\n';

            for (int k = 0; k <= 2*order-1; ++k) {
                interpTheta_lvl_th.push_back(
                    Interp::evalLagrangeBasis(theta, branchThetas, k));

                // if (level) std::cout << ith << ' ' << k << ' ' << interpTheta_lvl_th[k] << '\n';
            }

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
        std::vector<size_t> s_lvl;

        const int nph = phis[level].size();
        const int mph = phis[level+1].size();

        for (size_t iph = 0; iph < nph; ++iph) {
            realVec interpPhi_lvl_ph;
            const double phi = phis[level][iph];

            // Find idx of child phi nearest parent phi
            const int s = std::floor(mph * phi / (2.0 * PI));

            // Assemble child phis interpolating parent phi
            realVec branchPhis;
            for (int jph = s+1-order; jph <= s+order; ++jph)

                // Wrap jph if not in [0, mph-1]
                // size_t jph_wrapped = Math::wrapIdxToRange(jph, mph); 
                
                branchPhis.push_back(2.0*PI*jph/static_cast<double>(mph));

            for (size_t k = 0; k <= 2*order-1; ++k)
                interpPhi_lvl_ph.push_back(
                    Interp::evalTrigBasis(phi, branchPhis, k));
            
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

    constexpr double q = 3.5; // TODO: Optimize this

    normedDists = []() {
        realVec dists;

        for (int dz = 2; dz <= 3; ++dz)
            for (int dy = 0; dy <= 3; ++dy)
                for (int dx = 0; dx <= 3; ++dx)
                    dists.push_back(vec3d(dx, dy, dz).norm());

        std::sort(dists.begin(), dists.end());
        dists.erase(std::unique(dists.begin(), dists.end()), dists.end());

        return dists;
    } ();

    for (size_t level = 0; level <= maxLevel; ++level) {
        const auto L = Ls[level];
        const int nps = std::floor(q*L);
        const double nodeLeng = rootLeng / pow(2.0, level);

        std::vector<cmplxVec> transl_lvl;
        
        for (int iDist = 0; iDist < normedDists.size(); ++iDist) {
            const auto dr = normedDists[iDist];

            cmplxVec transl_lvl_k;

            for (int ips = 0; ips < nps; ++ips) {
                const double psi = PI * ips / (nps - 1.0); // or / nps ?
                cmplx coeff = 0.0;

                for (int l = 0; l <= L; ++l)
                    coeff +=
                        powI(l) * (2.0*l+1.0) *
                        sphericalHankel1(wavenum*dr*nodeLeng, l) *
                        legendreL(cos(psi), l).first;

                transl_lvl_k.push_back(coeff);
            }

            transl_lvl.push_back(transl_lvl_k);

        }

        transl.push_back(transl_lvl);
    }
};
