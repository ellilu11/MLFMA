#include "leaf.h"

LeafVec Leaf::leaves;

Leaf::Leaf(
    const SrcVec& srcs,
    const int branchIdx,
    Stem* const base)
    : Node(srcs, branchIdx, base)
{
    maxLevel = std::max(level, maxLevel);
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

    nodes.push_back(shared_from_this()); // 
}

/* buildRadPats()
 * Build radiation patterns due to sources in all leaves
 */
void Leaf::buildRadPats() {

    for (const auto& leaf : leaves) {

        const int level = leaf->level;
        const auto& center = leaf->center;

        const auto [nth, nph] = getNumAngles(level);

        for (int angIdx = 0; angIdx < nth*nph; ++angIdx) {

            const auto& kvec = tables.khat[level][angIdx] * wavenum;
            const auto& toThPh = tables.toThPh[level][angIdx];

            std::vector<vec2cd> radPat(leaf->srcs.size(), vec2cd::Zero());

            int srcIdx = 0;
            for (const auto& src : leaf->srcs)
                radPat[srcIdx++] = toThPh * src->getRadAlongDir(center, kvec);

            leaf->radPats.push_back(radPat);
        }
    }
}

/* buildMpoleCoeffs()
 * (S2M) Build multipole coefficients from sources in this node  
 */
void Leaf::buildMpoleCoeffs() {

    if (isSrcless() || isRoot()) return;

    const auto [nth, nph] = getNumAngles(level);

    auto start = Clock::now();

    for (int angIdx = 0; angIdx < nth*nph; ++angIdx) {

        vec2cd coeff = vec2cd::Zero();

        int srcIdx = 0;
        for (const auto& src : srcs)
            // coeff += src->getCurrent() * radPats[angIdx][srcIdx++];
            coeff += solver->getQvec(src->getIdx()) * radPats[angIdx][srcIdx++];

        coeffs.push_back(coeff);

    }

    t.S2M += Clock::now() - start;

}

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Leaf::buildLocalCoeffs() {
    if (isRoot()) return;

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

/* evalFarSols()
 * (L2T) Evaluate sols from local expansion due to far nodes
 */
void Leaf::evalFarSols() {
    if (isSrcless() || level <= 1) return;

    const auto [nth, nph] = getNumAngles(level);

    const double phiWeight = 2.0*PI / static_cast<double>(nph); // TODO: static member

    int obsIdx = 0;
    for (const auto& obs : srcs) {

        size_t angIdx = 0;
        cmplx sol = 0;

        for (int ith = 0; ith < nth; ++ith) {
            const double weight = thetaWeights[level][ith];

            for (int iph = 0; iph < nph; ++iph) {
                // Do the angular integration
                sol += weight 
                    * radPats[angIdx][obsIdx].dot(localCoeffs[angIdx]); // Hermitian dot!

                ++angIdx;
            }
        }

        // obs->addToSol(Phys::C * wavenum * phiWeight * sol);
        solver->addToSols(obs->getIdx(), Phys::C * wavenum * phiWeight * sol);

        ++obsIdx;
    }
}

/* evalNearNonNborSols()
 * (M2T/S2T) Evaluate sols from mpole expansion due to list 3 nodes
 */
void Leaf::evalNearNonNborSols() {
    // No reciprocity
    //for (const auto& node : nearNonNbors)
    //    evalPairSols(node);
    //return;

    // Do nothing! Contribution from list 3 node was 
    // already evaluated by Node::evalLeafIlistSols()
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
 * Sum solutions at all sources in all leaves 
 */ 
void Leaf::evaluateSols() {

    auto start = Clock::now();
    for (const auto& leaf : leaves)
        leaf->evalFarSols();
    t.L2T += Clock::now() - start;

    start = Clock::now();

    for (const auto& leaf : leaves)
        leaf->evalNearNonNborSols();

    for (const auto& pair : findNearNborPairs()) {
        auto [obsLeaf, srcLeaf] = pair;
        obsLeaf->evalPairSols(srcLeaf);
    }

    // No reciprocity
    //for (const auto& obsLeaf : leaves)
    //    for (const auto& srcLeaf : obsLeaf->nearNbors)
    //        obsLeaf->evalPairSols(srcLeaf);

    for (const auto& leaf : leaves)
        leaf->evalSelfSols();
    
    t.S2T += Clock::now() - start;
}
