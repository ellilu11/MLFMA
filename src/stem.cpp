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
  
        for (int jth = 0; jth < mth; ++jth) {
            for (int jph = 0; jph < mph; ++jph) {

                const auto& kvec = tables.kvec[level+1][l];

                shiftedCoeffs.push_back(
                    exp(iu*kvec.dot(shift)) * branchCoeffs[l++]);

            }
        }

        // Interpolate over theta
        std::vector<vec3cd> interpedCoeffs(nth*mph, vec3cd::Zero()); 

        size_t m = 0;
        for (int ith = 0; ith < nth; ++ith) {
            auto t = tables.ts[level][ith];

            for (int jph = 0; jph < mph; ++jph) {

                for (int jth = t+1-order, k = 0; jth <= t+order; ++jth, ++k) {

                    // Flip jph if not in [0, mth-1]
                    const int jth_flipped = Math::flipIdxToRange(jth, mth);

                    const bool outOfRange = jth != jth_flipped; // jth < 0 || jth >= mth;

                    int jph_shifted = jph;

                    // if theta \notin [0, pi] then if:
                    // phi \in (0, pi) add pi, phi \in (pi, 2pi) subtract pi
                    if (outOfRange)
                        jph_shifted += ((jph < mph/2) ? mph/2 : -mph/2);

                    const int m_shifted = jth_flipped*mph + jph_shifted;

                    interpedCoeffs[m] +=
                        tables.interpTheta[level][ith][k] * shiftedCoeffs[m_shifted];
                        // * Math::pm(outOfRange); // only for spherical components!

                    /*if (outOfRange && jph == mph-1)
                        std::cout
                            << ith << ' '
                            << thetas[level+1][jth_flipped] << ' '
                            << phis[level+1][jph_shifted] << ' '
                            << tables.interpTheta[level][ith][k] << '\n';*/

                }            
                
                m++;
            }
        }
        
        // Interpolate over phi
        size_t n = 0;
        for (int ith = 0; ith < nth; ++ith) { 

            for (int iph = 0; iph < nph; ++iph) { 
                auto s = tables.ss[level][iph]; // don't need to lookup for every ith

                for (int jph = s+1-order, k = 0; jph <= s+order; ++jph, ++k) {

                    // Wrap jph if not in [0, mph-1]
                    const int jph_wrapped = Math::wrapIdxToRange(jph, mph); 

                    coeffs[n] +=
                        tables.interpPhi[level][iph][k]
                        * interpedCoeffs[ith*mph + jph_wrapped];
                }

                n++;
            }
        }

        // Shift polar branch coeffs to center and merge (no interp)
        polarCoeffs.first += exp(iu*vec3d(0, 0, 1).dot(shift)) 
            * branch->getPolarCoeffs().first;

        polarCoeffs.second += exp(iu*vec3d(0, 0, -1).dot(shift))
            * branch->getPolarCoeffs().second;

    } // for (const auto& branch : branches)

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

