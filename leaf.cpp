#include "leaf.h"

LeafVec Leaf::leaves;

Leaf::Leaf(
    const RWGVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 * Also find all neighbor leaves of equal or lesser size (list 1)
 */
void Leaf::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);

        if (nbor != nullptr) {
            nbors.push_back(nbor);
            auto nbors = getNeighborsLeqSize(nbor, dir);
            nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
        }
    }

    assert(nbors.size() <= numDir);
}

/* buildLists()
 * Add self to list of leaves 
 * Find neighbor and interaction lists
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void Leaf::buildLists() {
    leaves.push_back(shared_from_this()); 

    if (isRoot()) return;
    
    buildNeighbors();

    buildInteractionList();

    pushSelfToNearNonNbors();
}

/* buildMpoleCoeffs()
 * (P2M) Build mpole expansions from particles in this node  
 */
void Leaf::buildMpoleCoeffs() {
}

/* propagateExpCoeffs()
 * (M2X) Convert mpole coeffs into exp coeffs
 * (X2X) Translate exp coeffs to nodes in all dirlists
 */
void Leaf::propagateExpCoeffs() {
}

/* buildLocalCoeffs()
 * (X2L) Receive incoming exp coeffs and add to local coeffs
 * (P2L) Add contribution from list 4 nodes t to local coeffs
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Leaf::buildLocalCoeffs() {

}

/* evalFarSols()
 * (L2P) Evaluate sols from local expansion due to far nodes
 */
void Leaf::evalFarSols() {
}

/* evalNearNonNborSols()
 * (M2P) Evaluate sols from mpole expansion due to list 3 nodes
 */
void Leaf::evalNearNonNborSols() {
}

/* findNearNborPairs()
 * From list of leaves, find all near neighbor leaf pairs
 */
std::vector<LeafPair> Leaf::findNearNborPairs(){
    std::vector<LeafPair> leafPairs;

    for (const auto& leaf : leaves) {
        for (const auto& nbor : leaf->nearNbors) {
            auto nborLeaf = dynamic_pointer_cast<Leaf>(nbor);
            if (leaf < nborLeaf)
                leafPairs.emplace_back(leaf, nborLeaf);
        }
    }

    return leafPairs;
}

/* evaluateSols()
 * Sum solutions at all particles in all leaves 
 */ 
void Leaf::evaluateSols() {
}