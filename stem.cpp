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

