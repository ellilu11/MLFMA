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
    const int nth = thetas[level].size();
    const int nph = phis[level].size();
    coeffs.resize(nth*nph, vec3cd::Zero());

    const int mth = thetas[level+1].size();
    const int mph = phis[level+1].size();

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();

        // Shift branch coeffs to center of this box
        const auto rvec = center - branch->getCenter();
        size_t l = 0;
        std::vector<vec3cd> shiftedBranchCoeffs;

        for (int jth = 0; jth < mth; ++jth){
            for (int jph = 0; jph < mph; ++jph) {

                const auto kvec = tables.kvec[level+1][l];
                shiftedBranchCoeffs.push_back(
                    Math::expI(kvec.dot(rvec)) * branchCoeffs[l++]
                );
            }
        }

        // Interpolate over theta
        std::vector<vec3cd> interpedBranchCoeffs(nth*mph, vec3cd::Zero()); 
        size_t m = 0;
        for (int ith = 0; ith < nth; ++ith) { // over parent theta to interpolate
            auto t = tables.T[level][ith];

            for (int jph = 0; jph < mph; ++jph) { // over child phi (uninterpolated)

                for (int jth = t+1-order; jth <= t+order; ++jth) { // over child theta interpolating parent theta
                    const size_t k = jth - (t+1-order);
                    interpedBranchCoeffs[m] +=
                        tables.interpTheta[level][ith][k]
                        * shiftedBranchCoeffs[jth*mph + jph];
                }

                m++;
            }
        }

        // Interpolate over phi
        size_t n = 0;
        for (int ith = 0; ith < nth; ++ith) { // over parent theta (interpolated)

            for (int iph = 0; iph < nph; ++iph) { // over parent phi to interpolate
                auto s = tables.S[level][ith]; // don't need to lookup for every ith

                for (int jph = s+1-order; jph <= s+order; ++jph) { // over child phi interpolating parent phi
                    const size_t k = jph - (s+1-order);
                    coeffs[n] +=
                        tables.interpPhi[level][iph][k]
                        * interpedBranchCoeffs[ith*mph + jph];
                }

                n++;
            }
        }
    }
}

/* propagateExpCoeffs()
 * (M2X) Convert mpole coeffs into outgoing exp coeffs
 * (X2X) Translate outgoing exp coeffs to nodes in all dirlists
 */ 
void Stem::propagateExpCoeffs() {
}

/* buildLocalCoeffs() 
 * (X2L) Receive incoming exp coeffs and add to local coeffs
 * (P2L) Add contribution from list 4 nodes t to local coeffs
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Stem::buildLocalCoeffs() {
   
}

