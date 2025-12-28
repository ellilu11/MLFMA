#include "leaf.h"

LeafVec Leaf::leaves;
std::vector<LeafPair> Leaf::nearPairs;

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

/* findNearNborPairs()
 * From list of leaves, find all near neighbor leaf pairs
 */
void Leaf::findNearNborPairs() {

    for (const auto& leaf : leaves) {

        for (const auto& nbor : leaf->nearNbors) {

            auto nborLeaf = dynamic_pointer_cast<Leaf>(nbor);
            if (leaf < nborLeaf)
                nearPairs.emplace_back(leaf, nborLeaf);

        }
    }
}

void Leaf::buildNearRads() {

    findNearNborPairs();

    // return;

    for (const auto& pair : nearPairs) {
        auto [obsLeaf, srcLeaf] = pair;
        assert(obsLeaf < srcLeaf);

        for (size_t obsIdx = 0; obsIdx < obsLeaf->srcs.size(); ++obsIdx) {
            for (size_t srcIdx = 0; srcIdx < srcLeaf->srcs.size(); ++srcIdx) {
                const auto obs = obsLeaf->srcs[obsIdx], src = srcLeaf->srcs[srcIdx];

                const cmplx rad = obs->getIntegratedRad(src);

                obsLeaf->nearRads.push_back(rad); 
                //srcLeaf->nearRads.push_back(rad);
            }
        }
        // std::cout << obsLeaf->nearRads.size() << '\n';
    }

    for (const auto& leaf : leaves) {
        for (size_t obsIdx = 1; obsIdx < leaf->srcs.size(); ++obsIdx) { // obsIdx = 0
            for (size_t srcIdx = 0; srcIdx < obsIdx; ++srcIdx) { // srcIdx <= obsIdx 
                const auto& obs = leaf->srcs[obsIdx], src = leaf->srcs[srcIdx];

                leaf->selfRads.push_back(obs->getIntegratedRad(src));

            }
        }
    }
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
        for (const auto& src : srcs) {
            // coeff += src->getCurrent() * radPats[angIdx][srcIdx++];
            // std::cout << getCurrent(src->getIdx()) << '\n';
            coeff += (*currents)[src->getIdx()] * radPats[angIdx][srcIdx++];
        }

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
        // addToSol(obs->getIdx(), Phys::C * wavenum * phiWeight * sol);
        (*sols)[obs->getIdx()] += Phys::C * wavenum * phiWeight * sol;

        ++obsIdx;
    }
}

/* evalNearNonNborSols()
 * (M2T/S2T) Evaluate sols from mpole expansion due to list 3 nodes
 */
void Leaf::evalNearNonNborSols() {
    // No reciprocity
    //for (const auto& node : nearNonNbors)
    //    evalNonNearPairSols(node);
    //return;

    // Do nothing! Contribution from list 3 node was 
    // already evaluated by Node::evalLeafIlistSols()
}

/* evalPairSols(srcNode)
 * (S2T) Evaluate sols at sources in this node due to sources in srcNode
 * and vice versa
 * srcNode : source node
 */
void Leaf::evalPairSols(const std::shared_ptr<Node> srcNode) { // TODO: Debug

    auto srcLeaf = dynamic_pointer_cast<Leaf>(srcNode);

    const int numObss = srcs.size(), numSrcs = srcLeaf->srcs.size();

    cmplxVec solAtObss(numObss, 0.0);
    cmplxVec solAtSrcs(numSrcs, 0.0);

    int pairIdx = 0;
    for (size_t obsIdx = 0; obsIdx < numObss; ++obsIdx) {
        for (size_t srcIdx = 0; srcIdx < numSrcs; ++srcIdx) {
            const auto obs = srcs[obsIdx], src = srcLeaf->srcs[srcIdx];

            const cmplx rad = nearRads[pairIdx++];

            solAtObss[obsIdx] += (*currents)[src->getIdx()] * rad;
            solAtSrcs[srcIdx] += (*currents)[obs->getIdx()] * rad;

        }
    }

    for (int n = 0; n < numObss; ++n)
        (*sols)[srcs[n]->getIdx()] += Phys::C * wavenum * solAtObss[n];

    for (int n = 0; n < numSrcs; ++n)
        (*sols)[srcLeaf->srcs[n]->getIdx()] += Phys::C * wavenum * solAtSrcs[n];
}

/* evalSelfSols()
 * (S2T) Evaluate sols at sources in this node due to other sources
 * in this node
 */
void Leaf::evalSelfSols() {

    const int numSrcs = srcs.size();

    cmplxVec solAtObss(numSrcs, 0.0);

    // TODO: Handle self-interactions
    int pairIdx = 0;
    for (size_t obsIdx = 1; obsIdx < numSrcs; ++obsIdx) { // obsIdx = 0
        for (size_t srcIdx = 0; srcIdx < obsIdx; ++srcIdx) { // srcIdx <= obsIdx 
            auto obs = srcs[obsIdx], src = srcs[srcIdx];

            const auto glObsIdx = obs->getIdx(), glSrcIdx = src->getIdx();

            const cmplx rad = selfRads[pairIdx++];

            solAtObss[obsIdx] += (*currents)[src->getIdx()] * rad;
            solAtObss[srcIdx] += (*currents)[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < numSrcs; ++n)
        (*sols)[srcs[n]->getIdx()] += Phys::C * wavenum * solAtObss[n];

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

    //for (const auto& leaf : leaves)
    //    leaf->evalNearNonNborSols();

    for (const auto& pair : nearPairs) {
        auto [obsLeaf, srcLeaf] = pair;
        obsLeaf->evalPairSolsDir(srcLeaf);
    }

    // No reciprocity
    //for (const auto& obsLeaf : leaves)
    //    for (const auto& srcLeaf : obsLeaf->nearNbors)
    //        obsLeaf->evalPairSols(srcLeaf);

    for (const auto& leaf : leaves)
        leaf->evalSelfSols();
    
    t.S2T += Clock::now() - start;
}