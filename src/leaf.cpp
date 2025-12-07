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

    coeffs.resize(nth*nph, vec2cd::Zero());

    /*auto triCoeff = [this](std::shared_ptr<Triangle> tri, mat3d ImKK, vec3d kvec, bool isPlus) {
        auto [nodes, weight] = tri->getQuads();
        auto coeff = vec3cd::Zero();

        for (const auto& quadNode : nodes)
            coeff += weight * ImKK * (rwg->getVplus() - quadNode)
            * Math::expI(kvec.dot(center-quadNode));
    }*/

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {
        for (int iph = 0; iph < nph; ++iph) {

            const auto ImKK = tables.ImKK[level][idx];
            const auto kvec = tables.kvec[level][idx];

            vec3cd dirCoeff = vec3cd::Zero();
            for (const auto& rwg : rwgs) { // TODO: Optimize this loop
                vec3cd rwgCoeff = vec3cd::Zero();

                auto triPlus = rwg->getTriPlus();
                auto vPlus = rwg->getVplus();
                auto [nodesPlus,weightPlus] = triPlus->getQuads();
                for (const auto& quadNode : nodesPlus)
                    rwgCoeff += weightPlus * Math::expI(kvec.dot(center-quadNode)) 
                                    * (vPlus - quadNode);
                
                auto triMinus = rwg->getTriMinus();
                auto vMinus = rwg->getVminus();
                auto [nodesMinus, weightMinus] = triMinus->getQuads();
                for (const auto& quadNode : nodesMinus)
                    rwgCoeff += weightMinus * Math::expI(kvec.dot(center-quadNode)) 
                                    * (quadNode - vMinus);
                
                dirCoeff += rwg->getCurrent() * rwg->getLeng() * rwgCoeff;
            }

            // Get theta and phi components
            coeffs[idx] = tables.matToThPh[level][idx] * (ImKK * dirCoeff);

            // std::cout << '(' << ith << ',' << iph << ") " << coeffs[idx] << '\n';

            idx++;
        }
    }
    // std::cout << '\n';
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
 * Sum solutions at all particles in all leaves 
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