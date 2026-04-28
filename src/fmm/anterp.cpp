#include "farfield.h"

/* addAnterpCoeffs(inCoeffs, outCoeffs, srcLvl, tgtLvl)
 * Add anterpolated coefficients from inCoeffs to outCoeffs
 * srcLvl and tgtLvl are levels of the source and target nodes, respectively
 * Anterpolation is done in two steps: first over phi, then over theta
 */
void FMM::Farfield::addAnterpCoeffs(
    const Coeffs& inCoeffs, Coeffs& outCoeffs, int srcLvl, int tgtLvl) const
{
    int order = config.interpOrder;

    auto [mth, mph] = levels[srcLvl].getNumAngles();
    auto [nth, nph] = levels[tgtLvl].getNumAngles();
    assert(!(nph%2)); // nph needs to be even

    const int tblLvl = std::min(srcLvl, tgtLvl);
    const auto& interpTheta = levels[tblLvl].interpTheta;
    const auto& interpPhi = levels[tblLvl].interpPhi;

    // Anterpolate over extended phi
    Coeffs innerCoeffs(mth*nph);
    for (int iph = 0; iph < mph; ++iph) {
        const auto [interp, nearIdx] = interpPhi[iph];

        for (int jph = -order; jph < nph+order; ++jph) {
            int k = jph - (nearIdx+1-order);

            // If iph \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            int jph_wrapped = Math::wrapIdxToRange(jph, nph);

            for (int ith = 0; ith < mth; ++ith) {
                size_t idxInner = ith*nph+jph_wrapped, idxIn = ith*mph+iph;

                innerCoeffs.theta[idxInner] += interp[k] * inCoeffs.theta[idxIn];
                innerCoeffs.phi[idxInner] += interp[k] * inCoeffs.phi[idxIn];
            }
        }
    }

    // Anterpolate over extended theta
    for (int ith = 0; ith < mth; ++ith) {
        const auto [interp, nearIdx] = interpTheta[ith];

        for (int jth = -order; jth < nth+order; ++jth) {
            int k = jth - (nearIdx+1-order);

            // If ith \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            int jth_flipped = Math::flipIdxToRange(jth, nth);

            bool outOfRange = jth != jth_flipped; // jth < 0 || jth >= nth;

            for (int jph = 0; jph < nph; ++jph) {
                int jph_shifted = jph;
                if (outOfRange) jph_shifted += ((jph < nph/2) ? nph/2 : -nph/2);

                size_t idxOut = jth_flipped*nph+jph_shifted, idxInner = ith*nph+jph;

                outCoeffs.theta[idxOut] +=
                    interp[k] * innerCoeffs.theta[idxInner] * Math::sign(outOfRange);
                outCoeffs.phi[idxOut] +=
                    interp[k] * innerCoeffs.phi[idxInner] * Math::sign(outOfRange);
            }
        }
    }
}