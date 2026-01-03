#include "stem.h"

Stem::Stem(
    const SrcVec& srcs,
    const int branchIdx,
    Stem* const base)
    : Node(srcs, branchIdx, base)
{
    // Assign every src to a branch based on src center relative to node center
    std::array<SrcVec, 8> branchSrcs;

    for (const auto& src : srcs)
        branchSrcs[Math::bools2Idx(src->getCenter() > center)].push_back(src);

    // Construct branch nodes
    for (size_t k = 0; k < branchSrcs.size(); ++k) {
        std::shared_ptr<Node> branch;

        if (branchSrcs[k].size() > config.maxNodeSrcs)
            branch = std::make_shared<Stem>(branchSrcs[k], k, this);
        else
            branch = std::make_shared<Leaf>(branchSrcs[k], k, this);

        branches.push_back(std::move(branch));
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

        if (nbor) nbors.push_back(nbor);
    }

    assert(nbors.size() <= numDir);
}

/* initNode()
 * Find neighbor and interaction lists.
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void Stem::initNode() {
    resizeCoeffs();

    if (!isRoot()) {
        buildNeighbors();

        buildInteractionList();

        pushSelfToNearNonNbors();
    }

    for (const auto& branch : branches)
        branch->initNode();
}

/* buildMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs 
 */
void Stem::buildMpoleCoeffs() {
    const int order = config.interpOrder;

    const auto [mth, mph] = getNumAngles(level+1);

    std::fill(coeffs.begin(), coeffs.end(), vec2cd::Zero());

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
        addInterpCoeffs(shiftedCoeffs, coeffs, level+1, level);

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
    addInterpCoeffs(shiftedCoeffs, outCoeffs, level, level+1);

    return outCoeffs;
}

void Stem::addInterpCoeffs(
    const std::vector<vec2cd>& inCoeffs, 
    std::vector<vec2cd>& outCoeffs, 
    int srcLvl, int tgtLvl)
{
    const int order = config.interpOrder;

    const auto [mth, mph] = getNumAngles(srcLvl);
    const auto [nth, nph] = getNumAngles(tgtLvl);

    assert(!(mph%2)); // mph needs to be even

    // Select which interp tables to use
    const auto& interpTheta = 
        (srcLvl > tgtLvl) ? tables.interpTheta : tables.invInterpTheta;
    const auto& interpPhi = 
        (srcLvl > tgtLvl) ? tables.interpPhi : tables.invInterpPhi;

    const int tblLvl = std::min(srcLvl, tgtLvl);

    // Interpolate over theta
    std::vector<vec2cd> innerCoeffs(nth*mph, vec2cd::Zero());

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
                    interp[k] * inCoeffs[m_shifted] * Math::sign(outOfRange); 
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
                    interp[k] * innerCoeffs[jthmph + iph_wrapped];
            }

            ++n;
        }
    }
}

/*void Stem::addAnterpCoeffs(
    const std::vector<vec3cd>& inCoeffs, std::vector<vec3cd>& outCoeffs, int level)
{
    const int order = config.interpOrder;

    const auto [nth, nph] = getNumAngles(level);
    const auto [mth, mph] = getNumAngles(level+1);

    assert(!(mph%2)); // mph needs to be even

    std::vector<vec3cd> outCoeffs(mth*mph, vec3cd::Zero());

    // Anterpolate over theta
    std::vector<vec3cd> innerCoeffs(mth*nph, vec3cd::Zero());

    size_t m = 0;
    for (int ith = 0; ith < mth; ++ith) { // over child thetas to anterpolate

        for (int jph = 0; jph < nph; ++jph) { // over parent phis (unanterpolated)

            for (int jth = 0; jth < nth; ++jth) { // over parent thetas anterpolating child thetas

                const int t = tables.ts[level][jth]; // TODO: don't need to lookup for every ith & jph

                // shift from ith \in [t+1-order,t+order] to k \in [0,2*order-1]   
                const int k = ith - (t+1-order);

                // std::cout << ith << ' ' << jth << ' ' << t << ' ' << k << '\n';

                // if ith \notin [t+1-order,t+order], matrix element is zero
                if (k < 0 || k >= 2*order) continue;

                innerCoeffs[m] +=
                    tables.interpTheta[level][jth][k] * inCoeffs[m];

                // std::cout << tables.interpTheta[level][jth][k] << ' ';
            }

            // std::cout << std::setprecision(9) << inCoeffs[m] << '\n';

            m++;
        }
    }

    // Anterpolate over phi
    size_t n = 0;
    for (int ith = 0; ith < mth; ++ith) { // over child thetas (anterpolated)

        for (int iph = 0; iph < mph; ++iph) { // over child phis to anterpolate

            for (int jph = 0; jph < nph; ++jph) { // over parent phis anterpolating child phis

                const int s = tables.ss[level][jph]; // TODO: don't need to lookup for every ith & iph

                // shift from iph \in [s+1-order,s+order] to k \in [0,2*order-1]
                const int k = iph - (s+1-order);

                // if iph \notin [s+1-order,s+order], matrix element is zero
                if (k < 0 || k >= 2*order) continue;

                outCoeffs[n] +=
                    tables.interpPhi[level][jph][k] * innerCoeffs[n];
            }

            // std::cout << std::setprecision(9) << outCoeffs[n] << '\n';

            n++;
        }
    }

    return outCoeffs;
}*/

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