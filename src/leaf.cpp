#include "leaf.h"

LeafVec Leaf::leaves;

Leaf::Leaf(
    const RWGVec& rwgs,
    const int branchIdx,
    Stem* const base)
    : Node(rwgs, branchIdx, base)
{
    maxLevel = level;
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
 * (S2M) Build multipole coefficients from RWG in this node  
 */
void Leaf::buildMpoleCoeffs() {
    if (rwgs.empty()) return;

    const auto [nth, nph] = getNumAngles(level);

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto& ImKK = tables.ImKK[level][idx];
            const auto& kvec = tables.kvec[level][idx];

            vec3cd dirCoeff = vec3cd::Zero();
            for (const auto& rwg : rwgs)
                dirCoeff += rwg->getRadAlongDir(center, kvec);

            // in spherical components
            // coeffs.push_back(tables.matToSph[level][idx] * (ImKK * dirCoeff));

            // in spherical (no radial) components
            // coeffs.push_back(tables.matToThPh[level][idx] * (ImKK * dirCoeff));

            // in cartesian components
            coeffs.push_back(ImKK * dirCoeff);

            idx++;
        }
    }

    // Get polar coeffs in cartesian components
    vec3cd northCoeff = vec3cd::Zero();
    for (const auto& rwg : rwgs)
        northCoeff += rwg->getRadAlongDir(center, vec3d(0,0,wavenum));
    polarCoeffs.first = Math::IminusRR(0, 0) * northCoeff;

    vec3cd southCoeff = vec3cd::Zero();
    for (const auto& rwg : rwgs)
        southCoeff += rwg->getRadAlongDir(center, vec3d(0,0,-wavenum));
    polarCoeffs.second = Math::IminusRR(0, PI) * southCoeff;

}

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Leaf::buildLocalCoeffs() {
    if (isRoot()) return;

    buildMpoleToLocalCoeffs();

    evalLeafIlistSols();

    if (!base->isRoot()) {
        auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);

        for (int l = 0; l <= order; ++l)
            localCoeffs[l] += shiftedLocalCoeffs[l];
    }
}

/* evalFarSols()
 * (L2T) Evaluate sols from local expansion due to far nodes
 */
void Leaf::evalFarSols() {
}

/* evalNearNonNborSols()
 * (M2T) Evaluate sols from mpole expansion due to list 3 nodes
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
 * Sum solutions at all RWGs in all leaves 
 */ 
void Leaf::evaluateSols() {

    for (const auto& leaf : leaves) {
        leaf->evalFarSols();

        leaf->evalNearNonNborSols();

        leaf->evalSelfSols();
    }

    for (const auto& pair : findNearNborPairs()) {
        auto [obsLeaf, srcLeaf] = pair;
        obsLeaf->evalPairSols(srcLeaf);
    }

}