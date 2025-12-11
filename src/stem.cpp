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
 * Find neighbor and interaction lists.,tr1qsawxdc2qw3eryuiol;qwsxc 
 * * Add self as near non-neighbor (list 3 node) of any list 4 nodes
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

    assert(!(mph%2)); // mph needs to be even

    coeffs.resize(nth*nph, vec3cd::Zero());

    for (const auto& branch : branches) {
        if (branch->getRWGs().empty()) continue;

        branch->buildMpoleCoeffs();
        const auto& branchCoeffs = branch->getMpoleCoeffs();

        // Shift branch coeffs to center of this box
        const auto& shift = center - branch->getCenter();

        std::vector<vec3cd> shiftedCoeffs;

        size_t l = 0;
        for (int ith = 0; ith < mth; ++ith) {
            for (int iph = 0; iph < mph; ++iph) {

                const auto& kvec = tables.kvec[level+1][l];

                shiftedCoeffs.push_back(
                    exp(iu*kvec.dot(shift)) * branchCoeffs[l++]);

            }
        }

        // Interpolate over theta
        std::vector<vec3cd> interpedCoeffs(nth*mph, vec3cd::Zero()); 

        size_t m = 0;
        for (int jth = 0; jth < nth; ++jth) {
            const int t = tables.ts[level][jth];

            for (int iph = 0; iph < mph; ++iph) {

                for (int ith = t+1-order, k = 0; ith <= t+order; ++ith, ++k) {

                    // Flip iph if not in [0, mth-1]
                    const int ith_flipped = Math::flipIdxToRange(ith, mth);

                    const bool outOfRange = ith != ith_flipped; // jth < 0 || jth >= mth;

                    int iph_shifted = iph;

                    // if theta \notin [0, pi] then if:
                    // phi \in (0, pi) add pi, phi \in (pi, 2pi) subtract pi
                    if (outOfRange)
                        iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                    const int m_shifted = ith_flipped*mph + iph_shifted;

                    interpedCoeffs[m] +=
                        tables.interpTheta[level][jth][k] * shiftedCoeffs[m_shifted];
                        // * Math::pm(outOfRange); // only for spherical components!
                }            
                
                m++;
            }
        }
        
        // Interpolate over phi
        size_t n = 0;
        for (int jth = 0; jth < nth; ++jth) { 

            for (int jph = 0; jph < nph; ++jph) { 
                const int s = tables.ss[level][jph]; // TODO: don't need to lookup for every ith

                for (int iph = s+1-order, k = 0; iph <= s+order; ++iph, ++k) {

                    // Wrap iph if not in [0, mph-1]
                    const int iph_wrapped = Math::wrapIdxToRange(iph, mph); 

                    coeffs[n] +=
                        tables.interpPhi[level][jph][k]
                        * interpedCoeffs[jth*mph + iph_wrapped];
                }

                n++;
            }
        }

        // Shift polar branch coeffs to center and merge (no interp)
        polarCoeffs.first += exp(iu*northPole.dot(shift)) 
            * branch->getPolarCoeffs().first;

        polarCoeffs.second += exp(iu*southPole.dot(shift))
            * branch->getPolarCoeffs().second;

    } // for (const auto& branch : branches)

    /*for (int ith = 0; ith < nth; ++ith)
        for (int iph = 0; iph < nph; ++iph)
            std::cout << '(' << level << ',' << ith << ',' << iph << ") " << coeffs[ith*nph+iph] << '\n';*/
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
std::vector<vec3cd> Stem::getShiftedLocalCoeffs(const int branchIdx) const {

    const auto [nth, nph] = getNumAngles(level);
    const auto [mth, mph] = getNumAngles(level+1);

    std::vector<vec3cd> outCoeffs(mth*mph, vec3cd::Zero());

    // Shift local coeffs to center of branch
    const auto& shift = branches[branchIdx]->getCenter() - center;

    std::vector<vec3cd> shiftedCoeffs;

    size_t l = 0;
    for (int jth = 0; jth < nth; ++jth) {
        for (int jph = 0; jph < nph; ++jph) {

            const auto& kvec = tables.kvec[level][l];

            shiftedCoeffs.push_back(
                exp(iu*kvec.dot(shift)) * localCoeffs[l++]);

        }
    }

    // Anterpolate over theta
    std::vector<vec3cd> anterpedCoeffs(mth*nph, vec3cd::Zero());

    size_t m = 0;
    for (int ith = 0; ith < mth; ++ith) { // over child thetas to anterpolate

        for (int jph = 0; jph < nph; ++jph) { // over parent phis (unanterpolated)

            for (int jth = 0; jth < nth; ++jth) { // over parent thetas anterpolating child thetas

                const int t = tables.ts[level][jth]; // TODO: don't need to lookup for every ith & jph
                
                // shift from ith \in [t+1-order,t+order] to k \in [0,2*order-1]   
                const int k = ith - (t+1-order); 
                
                if (k < 0 || k >= 2*order) continue;
  
                anterpedCoeffs[m] +=
                    tables.interpTheta[level][jth][k] * shiftedCoeffs[m];
            }

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

                if (k < 0 || k >= 2*order) continue;

                outCoeffs[n] +=
                    tables.interpPhi[level][jph][k] * anterpedCoeffs[n];
            }

            n++;
        }
    }

    return outCoeffs;
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
            auto baseStem = static_cast<Stem*>(base);
                
            auto shiftedLocalCoeffs = baseStem->getShiftedLocalCoeffs(branchIdx);

            for (int l = 0; l <= order; ++l)
                localCoeffs[l] += shiftedLocalCoeffs[l];
        }
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

