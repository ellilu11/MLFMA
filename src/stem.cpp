#include "stem.h"

Stem::Stem(
    const RWGVec& rwgs,
    const int branchIdx,
    Stem* const base)
    : Node(rwgs, branchIdx, base)
{
    // Assign every RWG to a branch based on RWG center relative to node center
    std::array<RWGVec,8> branchRWGs;

    for (const auto& rwg : rwgs)
        branchRWGs[Math::bools2Idx(rwg->getCenter() > center)].push_back(rwg);
 
    // Construct branch nodes
    for (size_t k = 0; k < branchRWGs.size(); ++k) {
        std::shared_ptr<Node> branch;

        if (branchRWGs[k].size() > maxNodeSrcs)
            branch = std::make_shared<Stem>(branchRWGs[k], k, this);
        else
            branch = std::make_shared<Leaf>(branchRWGs[k], k, this);

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
        if (nbor != nullptr)
            nbors.push_back(nbor);
    }
    assert(nbors.size() <= numDir);
}

/* buildLists()
 * Find neighbor and interaction lists
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void Stem::buildLists() {
    if (!isRoot()) {
        buildNeighbors();

        buildInteractionList();

        pushSelfToNearNonNbors();
    }

    for (const auto& branch : branches)
        branch->buildLists();
}

/* buildMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs 
 */
void Stem::buildMpoleCoeffs() {
    const auto [nth, nph] = getNumAngles(level);
    const auto [mth, mph] = getNumAngles(level+1);

    coeffs.resize(nth*nph, vec2cd::Zero());

    for (const auto& branch : branches) {
        if (branch->getRWGs().empty()) continue;

        branch->buildMpoleCoeffs();

        const auto branchCoeffs = branch->getMpoleCoeffs();

        // Shift branch coeffs to center of this box
        const auto shift = center - branch->getCenter();

        std::vector<vec2cd> shiftedBranchCoeffs;

        size_t l = 0;
        for (int jth = 0; jth < mth; ++jth){
            for (int jph = 0; jph < mph; ++jph) {

                const auto kvec = tables.kvec[level+1][l];
                shiftedBranchCoeffs.push_back(
                    Math::expI(kvec.dot(shift)) * branchCoeffs[l++]);
            }
        }

        // Interpolate over theta
        std::vector<vec2cd> interpedBranchCoeffs(nth*mph, vec2cd::Zero()); 

        size_t m = 0;
        for (int ith = 0; ith < nth; ++ith) {
            auto t = tables.ts[level][ith];

            for (int jph = 0; jph < mph; ++jph) {

                size_t k = 0;
                for (int jth = t+1-order; jth <= t+order; ++jth) { 
                    // const int k = jth - (t+1-order);

                    // Flip jph if not in [0, mth-1]
                    const int jth_flipped = Math::flipIdxToRange(jth, mth);

                    bool outOfRange = jth != jth_flipped;// jth < 0 || jth >= mth;

                    int m2 = jth_flipped*mph + jph;
                    if (outOfRange)
                        m2 += (jph < mph/2) ? mph/2 : -mph/2;

                    interpedBranchCoeffs[m] +=
                        tables.interpTheta[level][ith][k++]
                        * shiftedBranchCoeffs[m2]
                        * Math::pm(outOfRange);
                }   

                // if (level)
                //    std::cout << '(' << ith << ',' << jph << ") " << interpedBranchCoeffs[m] << '\n';

                m++;
            }
        }
        
        // Interpolate over phi
        size_t n = 0;
        for (int ith = 0; ith < nth; ++ith) { 

            for (int iph = 0; iph < nph; ++iph) { 
                auto s = tables.ss[level][iph]; // don't need to lookup for every ith

                size_t k = 0;
                for (int jph = s+1-order; jph <= s+order; ++jph) {
                    // const int k = jph - (s+1-order);

                    // Wrap jph if not in [0, mph-1]
                    const int jph_wrapped = Math::wrapIdxToRange(jph, mph); 

                    coeffs[n] +=
                        tables.interpPhi[level][iph][k++]
                        * interpedBranchCoeffs[ith*mph + jph_wrapped];
                }

                n++;
            }
        }
    }

    /*for (int ith = 0; ith < nth; ++ith)
        for (int iph = 0; iph < nph; ++iph)
            std::cout << '(' << level << ',' << ith << ',' << iph << ") " << coeffs[ith*nph+iph] << '\n';*/
}

/* buildLocalCoeffs() 
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Stem::buildLocalCoeffs() {
    if (!isRoot()) {

        buildMpoleToLocalCoeffs();

        evalLeafIlistSols();

        if (!base->isRoot()) {
            auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);

            for (int l = 0; l <= order; ++l)
                localCoeffs[l] += shiftedLocalCoeffs[l];
        }
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

