#include "farfield.h"

/* addInterpCoeffs(inCoeffs, outCoeffs, srcLvl, tgtLvl)
 * Add interpolated coefficients from inCoeffs to outCoeffs
 * srcLvl and tgtLvl are levels of the source and target nodes, respectively
 * Interpolation is done in two steps: first over theta, then over phi
 */
void FMM::Farfield::addInterpCoeffs(
    const Coeffs& inCoeffs, Coeffs& outCoeffs, int srcLvl, int tgtLvl) const
{
    int order = config.interpOrder;

    auto [mth, mph] = levels[srcLvl].getNumAngles();
    auto [nth, nph] = levels[tgtLvl].getNumAngles();
    assert(!(mph%2)); // mph needs to be even

    int tblLvl = std::min(srcLvl, tgtLvl);
    const auto& interpTheta = levels[tblLvl].interpTheta;
    const auto& interpPhi = levels[tblLvl].interpPhi;

    // Interpolate over theta
    Coeffs innerCoeffs(nth*mph);
    for (int jth = 0; jth < nth; ++jth) {
        const auto [interp, nearIdx] = interpTheta[jth];

        for (int ith = nearIdx+1-order; ith <= nearIdx+order; ++ith) {
            int ith_flipped = Math::flipIdxToRange(ith, mth);
            int k = ith - (nearIdx+1-order);

            bool outOfRange = ith != ith_flipped; // jth < 0 || jth >= mth;

            for (int iph = 0; iph < mph; ++iph) {
                int iph_shifted = iph;
                if (outOfRange) iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                size_t idxInner = jth*mph+iph, idxIn = ith_flipped*mph+iph_shifted;

                innerCoeffs.theta[idxInner] +=
                    interp[k] * inCoeffs.theta[idxIn] * Math::sign(outOfRange);
                innerCoeffs.phi[idxInner] +=
                    interp[k] * inCoeffs.phi[idxIn] * Math::sign(outOfRange);
            }
        }
    }

    // Interpolate over phi
    for (int jph = 0; jph < nph; ++jph) {
        const auto [interp, nearIdx] = interpPhi[jph];

        for (int iph = nearIdx+1-order; iph <= nearIdx+order; ++iph) {
            int iph_wrapped = Math::wrapIdxToRange(iph, mph);
            int k = iph - (nearIdx+1-order);

            for (int jth = 0; jth < nth; ++jth) {
                size_t idxOut = jth*nph+jph, idxInner = jth*mph+iph_wrapped;

                outCoeffs.theta[idxOut] += interp[k] * innerCoeffs.theta[idxInner];
                outCoeffs.phi[idxOut] += interp[k] * innerCoeffs.phi[idxInner];
            }
        }
    }
}

/* addInterpCoeffAlongTh()
 * Add interpolated coefficients along jth to all nodes in nodes, 
 * from srcLvl to tgtLvl
 */
void FMM::Farfield::addInterpCoeffAlongTh(
    const NodeVec& nodes, int jth, int srcLvl, int tgtLvl) const 
{
    int order = config.interpOrder;

    auto [mth, mph] = levels[srcLvl].getNumAngles();
    auto [nth, nph] = levels[tgtLvl].getNumAngles();
    size_t mDir = mth*mph;
    assert(!(mph%2)); // mph needs to be even

    int tblLvl = std::min(srcLvl, tgtLvl);
    const auto& interpTheta = levels[tblLvl].interpTheta;
    const auto& interpPhi = levels[tblLvl].interpPhi;

    for (auto& node : nodes) {
        for (const auto& branch : node->branches) {
            if (branch->isSrcless()) continue;

            // Get shifted coeffs of this branch
            vec3d dX = node->center - branch->center;

            Coeffs shiftedCoeffs(mDir);
            for (int iDir = 0; iDir < mDir; ++iDir) {
                const vec3d& kvec = levels[srcLvl].khat[iDir] * config.k;
                cmplx shift = exp(iu*kvec.dot(dX));

                shiftedCoeffs.theta[iDir] = shift * branch->coeffs.theta[iDir];
                shiftedCoeffs.phi[iDir] = shift * branch->coeffs.phi[iDir];
            }

            // Get interpolation coefficients and indices of angles nearest target theta
            const auto [interpTh, nearTh] = interpTheta[jth];

            // Assemble source thetas and phis interpolating target theta and phi
            Coeffs innerCoeffs(2*order*nph);

            // Interpolate over phi
            for (int jph = 0; jph < nph; ++jph) {
                const auto [interpPh, nearPh] = interpPhi[jph];

                for (int iph = nearPh+1-order; iph <= nearPh+order; ++iph) {
                    int iph_wrapped = Math::wrapIdxToRange(iph, mph);
                    int k = iph - (nearPh+1-order);

                    for (int ith = nearTh+1-order; ith <= nearTh+order; ++ith) {
                        size_t idxIn = ith*mph+iph_wrapped;
                        size_t idxInner = (ith-(nearTh+1-order))*mph+jph;

                        innerCoeffs.theta[idxInner] += interpPh[k] * shiftedCoeffs.theta[idxIn];
                        innerCoeffs.phi[idxInner] += interpPh[k] * shiftedCoeffs.phi[idxIn];
                    }
                }
            }

            // Interpolate over theta
            for (int ith = nearTh+1-order; ith <= nearTh+order; ++ith) {
                int ith_flipped = Math::flipIdxToRange(ith, mth);
                int k = ith - (nearTh+1-order);
                bool outOfRange = ith != ith_flipped; // jth < 0 || jth >= mth;

                for (int jph = 0; jph < nph; ++jph) {
                    int jph_shifted = jph;
                    if (outOfRange) jph_shifted += ((jph < nph/2) ? nph/2 : -nph/2);

                    size_t idxInner = (ith-(nearTh+1-order))*mph+jph_shifted;
                    size_t idxOut = jth*nph+jph;

                    node->coeffs.theta[idxOut] +=
                        interpTh[k] * innerCoeffs.theta[idxInner] * Math::sign(outOfRange);
                    node->coeffs.phi[idxOut] +=
                        interpTh[k] * innerCoeffs.phi[idxInner] * Math::sign(outOfRange);
                }
            }
        }
    }
}