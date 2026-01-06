#include "stem.h"

Stem::Stem(
    const SrcVec& srcs,
    const int branchIdx,
    Stem* const base)
    : Node(srcs, branchIdx, base)
{
    // Assign every src to a branch based on src center relative to node center
    std::array<SrcVec,8> branchSrcs;

    for (const auto& src : srcs)
        branchSrcs[Math::bools2Idx(src->getCenter() > center)].push_back(src);
 
    // Construct branch nodes
    for (size_t k = 0; k < branchSrcs.size(); ++k) {
        std::shared_ptr<Node> branch;

        if (branchSrcs[k].size() > config.maxNodeSrcs)
            branch = std::make_shared<Stem>(branchSrcs[k], k, this);
        else
            branch = std::make_shared<Leaf>(branchSrcs[k], k, this);

        branches.push_back(branch);
    }
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 */
void Stem::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);

        if (nbor != nullptr) nbors.push_back(nbor);
    }

    assert(nbors.size() <= numDir);
}

/* buildLists()
 * Find neighbor and interaction lists.
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void Stem::buildLists() {
    if (!isRoot()) {
        buildNeighbors();

        buildInteractionList();

        pushSelfToNearNonNbors();

        nodes.push_back(shared_from_this()); // 
    }

    for (const auto& branch : branches)
        branch->buildLists();
}

/* buildMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs
 */
void Stem::buildMpoleCoeffs() {

    const int order = config.interpOrder;

    const auto [mth, mph] = getNumAngles(level+1);
    const auto [nth, nph] = getNumAngles(level);

    coeffs.resize(nth*nph, vec2cd::Zero());

    for (const auto& branch : branches) {
        if (branch->isSrcless()) continue;

        branch->buildMpoleCoeffs();

        auto start = Clock::now();

        const auto& branchCoeffs = branch->getMpoleCoeffs();

        // Shift branch coeffs to center of this node
        const auto& dX = center - branch->getCenter();

        std::vector<vec2cd> shiftedCoeffs(mth*mph, vec2cd::Zero());

        for (int idx = 0; idx < mth*mph; ++idx) {
            const auto& kvec = tables.khat[level+1][idx] * wavenum;

            shiftedCoeffs[idx] = exp(iu*kvec.dot(dX)) * branchCoeffs[idx];
        }

        // Interpolate shifted coeffs to this node's angular grid
        addInterpCoeffs<vec2cd>(shiftedCoeffs, coeffs, level+1, level);

        t.M2M += Clock::now() - start;

    }
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
std::vector<vec2cd> Stem::getShiftedLocalCoeffs(int branchIdx) const {

    const auto [mth, mph] = getNumAngles(level);
    const auto [nth, nph] = getNumAngles(level+1);

    std::vector<vec2cd> outCoeffs(nth*nph, vec2cd::Zero());
    if (iList.empty()) return outCoeffs;

    // Shift local coeffs to center of branch
    const auto& dX = branches[branchIdx]->getCenter() - center;

    std::vector<vec2cd> shiftedCoeffs(mth*mph, vec2cd::Zero());

    for (int idx = 0; idx < mth*mph; ++idx) {
        const auto& kvec = tables.khat[level][idx] * wavenum;

        shiftedCoeffs[idx] = exp(iu*kvec.dot(dX)) * localCoeffs[idx];
    }

    // Interpolate shifted coeffs to branch's angular grid
    addAnterpCoeffs<vec2cd>(shiftedCoeffs, outCoeffs, level, level+1);

    return outCoeffs;
}

template <typename T>
void Stem::addInterpCoeffs(
    const std::vector<T>& inCoeffs, std::vector<T>& outCoeffs, int srcLvl, int tgtLvl)
{
    const int order = config.interpOrder;

    const auto [mth, mph] = getNumAngles(srcLvl);
    const auto [nth, nph] = getNumAngles(tgtLvl);

    assert(!(mph%2)); // mph needs to be even

    // Choose which interp tables to use
    const auto& interpTheta =
        (srcLvl > tgtLvl) ? tables.interpTheta : tables.invInterpTheta;
    const auto& interpPhi =
        (srcLvl > tgtLvl) ? tables.interpPhi : tables.invInterpPhi;

    const int tblLvl = std::min(srcLvl, tgtLvl);

    // Interpolate over theta
    std::vector<T> innerCoeffs(nth*mph, T{});

    size_t m = 0;
    for (int jth = 0; jth < nth; ++jth) {

        const auto [interp, nearIdx] = interpTheta[tblLvl][jth];

        for (int iph = 0; iph < mph; ++iph) {

            for (int ith = nearIdx+1-order, k = 0; ith <= nearIdx+order; ++ith, ++k) {

                // Flip ith if not in [0, mth-1]
                const int ith_flipped = Math::flipIdxToRange(ith, mth);

                const bool outOfRange = ith != ith_flipped; // jth < 0 || jth >= mth;

                int iph_shifted = iph;

                // if theta \notin [0, pi] then if:
                // phi \in (0, pi) add pi, phi \in (pi, 2pi) subtract pi
                if (outOfRange)
                    iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                const int m_shifted = ith_flipped*mph + iph_shifted;

                innerCoeffs[m] +=
                    interp[k] * inCoeffs[m_shifted]
                    * Math::sign(outOfRange); // only for spherical vector components!
            }

            ++m;
        }
    }

    // Interpolate over phi
    size_t n = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const int jthmph = jth*mph;

        for (int jph = 0; jph < nph; ++jph) {
            const auto [interp, nearIdx] = interpPhi[tblLvl][jph]; // TODO: don't need to lookup for every jth

            for (int iph = nearIdx+1-order, k = 0; iph <= nearIdx+order; ++iph, ++k) {

                // Wrap iph if not in [0, mph-1]
                const int iph_wrapped = Math::wrapIdxToRange(iph, mph);

                outCoeffs[n] +=
                    interp[k]
                    * innerCoeffs[jthmph + iph_wrapped];
            }

            ++n;
        }
    }
}

template <typename T>
void Stem::addAnterpCoeffs(
    const std::vector<T>& inCoeffs, std::vector<T>& outCoeffs, int srcLvl, int tgtLvl)
{
    const int order = config.interpOrder;

    const auto [mth, mph] = getNumAngles(srcLvl);
    const auto [nth, nph] = getNumAngles(tgtLvl);
    const int Nth = nth+2*order, Nph = nph+2*order;

    assert(!(nph%2)); // nph needs to be even

    const int tblLvl = std::min(srcLvl, tgtLvl);

    // Anterpolate over phi
    std::vector<T> innerCoeffs(mth*Nph, T{});
    for (int iph = 0; iph < mph; ++iph) { // over parent phis anterpolating child phis
        const auto [interp, nearIdx] = tables.interpPhi[tblLvl][iph];

        for (int jph = -order; jph < nph+order; ++jph) { // over child phis to anterpolate

            // shift from jph \in [nearIdx+1-order,nearIdx+order] to k \in [0,2*order-1]
            const int k = jph - (nearIdx+1-order);

            // if iph \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            for (int ith = 0; ith < mth; ++ith) // over parent thetas (unanterpolated)
                innerCoeffs[ith*Nph+jph+order] += interp[k] * inCoeffs[ith*mph+iph];
        }
    }

    // Anterpolate over theta
    std::vector<T> longCoeffs(Nth*Nph, T{});
    for (int ith = 0; ith < mth; ++ith) { // over parent thetas anterpolating child thetas

        const auto [interp, nearIdx] = tables.interpTheta[tblLvl][ith];
        for (int jth = -order; jth < nth+order; ++jth) { // over child thetas to anterpolate

            // shift from ith \in [t+1-order,t+order] to k \in [0,2*order-1]   
            const int k = jth - (nearIdx+1-order);

            // if ith \notin [t+1-order,t+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            for (int jph = -order; jph < nph+order; ++jph) // over child phis (anterpolated)
                longCoeffs[(jth+order)*Nph+jph+order] += interp[k] * innerCoeffs[ith*Nph+jph+order];
        }
    }

    // Contract grid points in theta and phi
    int dirIdx = 0;

    for (int jth = 0; jth < nth; ++jth) {
        for (int jph = 0; jph < nph; ++jph) {
            const int Jth = jth+order, Jph = jph+order;

            outCoeffs[dirIdx] = longCoeffs[Jth*Nph+Jph];

            // Handle nodes near prime meridian
            if (jph < order || jph >= nph-order) {
                const int Jph_wrapped = Jph + (jph < order ? nph : -nph);
                outCoeffs[dirIdx] += longCoeffs[(jth+order)*Nph+Jph_wrapped];
            }

            // Handle nodes near poles
            // TODO: Fix
            if (jth < order || jth >= nth-order) {
                const int Jph_shifted = Jph + ((jph < nph/2) ? nph/2 : -nph/2);
                const int jth_flipped = (jth < order ? -jth-1 : 2*nth-jth-1);

                outCoeffs[dirIdx] += longCoeffs[(jth_flipped+order)*Nph+Jph_shifted]
                    * Math::sign(1); // spherical vector components only
            }

            ++dirIdx;
        }
    }
}
//

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Stem::buildLocalCoeffs() {
    if (!isRoot()) {

        auto start = Clock::now();
        buildMpoleToLocalCoeffs();
        t.M2L += Clock::now() - start;

        evalLeafIlistSols();

        start = Clock::now();
        if (!base->isRoot()) {
            auto stemBase = dynamic_cast<Stem*>(base);

            localCoeffs =
                localCoeffs + stemBase->getShiftedLocalCoeffs(branchIdx);
        }

        t.L2L += Clock::now() - start;
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}