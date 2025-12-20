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

/* buildMpoleCoeffs()
 * (S2M) Build multipole coefficients from sources in this node  
 */
void Leaf::buildMpoleCoeffs() {

    if (isSrcless() || isRoot()) return;

    const auto [nth, nph] = getNumAngles(level);

    auto start = Clock::now();

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto& kvec = tables.khat[level][idx] * wavenum;
            const auto& ImKK = tables.ImKK[level][idx];

            vec3cd radPat = vec3cd::Zero();

            for (const auto& src : srcs)
                radPat += src->getRadAlongDir(center, kvec);

            // in spherical components
            // coeffs.push_back(tables.matToSph[level][idx] * (ImKK * dirCoeff));

            // in spherical (no radial) components
            // coeffs.push_back(tables.matToThPh[level][idx] * (ImKK * dirCoeff));

            // in cartesian components
            coeffs.push_back(ImKK * radPat);

            idx++;
        }
    }

    /* Get polar coeffs in cartesian components
    vec3cd northCoeff = vec3cd::Zero();
    for (const auto& rwg : rwgs)
        northCoeff += rwg->getRadAlongDir(center, wavenum*northPole);
    polarCoeffs.first = Math::IminusRR(0, 0) * northCoeff;

    vec3cd southCoeff = vec3cd::Zero();
    for (const auto& rwg : rwgs)
        southCoeff += rwg->getRadAlongDir(center, wavenum*southPole);
    polarCoeffs.second = Math::IminusRR(0, PI) * southCoeff;
    */

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
        auto stemBase = static_cast<Stem*>(base);

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

    for (const auto& obs : srcs) {

        size_t idx = 0;

        cmplx sol = 0;

        for (int ith = 0; ith < nth; ++ith) {

            const auto theta = thetas[level][ith];

            for (int iph = 0; iph < nph; ++iph) {

                const auto& khat = tables.khat[level][idx];

                // Compute incoming pattern along khat at this source
                const auto& incPat =
                    // tables.ImKK[level][idx] * // Tangential-T
                    // -khat.cross( // Tangential-K
                    obs->getIncAlongDir(center, wavenum*khat);

                // Do the angular integration
                sol += thetaWeights[level][ith]
                        * sin(theta) // TODO: Absorb into thetaWeights
                        * incPat.dot(localCoeffs[idx]);

                idx++;
            }
        }

        obs->addToSol(C * wavenum * phiWeight * sol);
    }

}

/* evalNearNonNborSols()
 * (M2T) Evaluate sols from mpole expansion due to list 3 nodes
 */
void Leaf::evalNearNonNborSols() {
    for (const auto& node : nearNonNbors)
        evalPairSols(node);
    return;
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

    auto start = Clock::now();

    for (const auto& leaf : leaves)
        leaf->evalFarSols();

    t.L2T += Clock::now() - start;

    for (const auto& leaf : leaves)
        leaf->evalNearNonNborSols();

    start = Clock::now();

    // Using reciprocity
    //for (const auto& pair : findNearNborPairs()) {
    //    auto [obsLeaf, srcLeaf] = pair;
    //    obsLeaf->evalPairSols(srcLeaf);
    //}

    for (const auto& obsLeaf : leaves)
        for (const auto& srcLeaf : obsLeaf->nearNbors)
            obsLeaf->evalPairSols(srcLeaf);

    for (const auto& leaf : leaves)
        leaf->evalSelfSols();

    t.S2T += Clock::now() - start;
}
