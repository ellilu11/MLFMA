#include "leaf.h"

LeafVec FMM::Leaf::leaves;
std::vector<LeafPair> FMM::Leaf::nearPairs;

FMM::Leaf::Leaf(
    const SrcVec& srcs,
    const int branchIdx,
    Stem* const base)
    : Node(srcs, branchIdx, base), 
      leafPairIdx(0), nonNearPairIdx(0)
{
    maxLevel = std::max(level, maxLevel);

    /* Assign indices to all sources in this leaf
    for (const auto& src : srcs) {
        src->setIdx(glSrcIdx++);
        std::cout << src->getIdx() << ' ';
    }
    if (!isSrcless()) std::cout << '\n';
    */
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 * Also find all neighbor leaves of equal or lesser size (list 1)
 */
void FMM::Leaf::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);

        if (nbor) {
            nbors.push_back(nbor);
            auto nbors = getNeighborsLeqSize(nbor, dir);
            nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
        }
    }

    assert(nbors.size() <= numDir);
}

/* initNode()
 * Add self to list of leaves 
 * Find neighbor and interaction lists
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void FMM::Leaf::initNode() {
    leaves.push_back(shared_from_this()); 

    resizeCoeffs();

    if (isRoot()) return;
    
    buildNeighbors();

    buildInteractionList();

    pushSelfToNearNonNbors();
}

/* findNearNborPairs()
 * From list of leaves, find all near neighbor leaf pairs
 */
void FMM::Leaf::findNearNborPairs() {

    for (const auto& leaf : leaves) {
        for (const auto& nbor : leaf->nearNbors) {
            auto nborLeaf = dynamic_pointer_cast<Leaf>(nbor);

            if (leaf < nborLeaf)
                nearPairs.emplace_back(leaf, nborLeaf);
        }
    }
}

void FMM::Leaf::buildNearRads() {

    findNearNborPairs();

    for (const auto& [obsLeaf, srcLeaf] : nearPairs) {
        assert(obsLeaf < srcLeaf);

        const size_t nObs = obsLeaf->srcs.size(), nSrcs = srcLeaf->srcs.size();

        auto leafPairRads = cmplxVec(nObs*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObs; ++iObs) {
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                const auto obs = obsLeaf->srcs[iObs], src = srcLeaf->srcs[iSrc];

                leafPairRads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }

        obsLeaf->nearRads.push_back(leafPairRads);
    }

    for (const auto& [obsNode, srcNode] : nonNearPairs) {
        auto obsLeaf = dynamic_pointer_cast<Leaf>(obsNode);

        const size_t nObs = obsLeaf->srcs.size(), nSrcs = srcNode->getSrcs().size();

        auto nodePairRads = cmplxVec(nObs*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObs; ++iObs) {
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                const auto obs = obsLeaf->srcs[iObs], src = srcNode->getSrcs()[iSrc];

                nodePairRads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }

        obsLeaf->nonNearRads.push_back(nodePairRads);
    }

    //std::ofstream zmatFile("out/zmat.txt");
    //std::ofstream vvecFile("out/vvec.txt");

    for (const auto& leaf : leaves) {
        for (size_t iObs = 1; iObs < leaf->srcs.size(); ++iObs) { // iObs = 0
            const auto& obs = leaf->srcs[iObs];

            for (size_t iSrc = 0; iSrc < iObs; ++iSrc) { // iSrc <= iObs 
                const auto& src = leaf->srcs[iSrc];

                leaf->selfRads.push_back(obs->getIntegratedRad(src));

            }
        }

        /* GMRES testing
        for (size_t iObs = 0; iObs < leaf->srcs.size(); ++iObs) {
            const auto& obs = leaf->srcs[iObs];

            for (size_t iSrc = 0; iSrc < leaf->srcs.size(); ++iSrc) {
                const auto& src = leaf->srcs[iSrc];

                zmatFile << Phys::C * wavenum * obs->getIntegratedRad(src) << ' ';

            }
            vvecFile << obs->getVoltage() << '\n';
            zmatFile << '\n';
        }
        */
    }
}

/* buildRadPats()
 * Build radiation patterns due to sources in all leaves
 */
void FMM::Leaf::buildRadPats() {

    for (const auto& leaf : leaves) {
        const int level = leaf->level;
        const auto& center = leaf->center;

        const auto [nth, nph] = angles[level].getNumAngles();
        const auto& tables_lvl = tables[level];

        for (int iDir = 0; iDir < nth*nph; ++iDir) {
            const auto& kvec = tables_lvl.khat[iDir] * wavenum;
            const auto& toThPh = tables_lvl.toThPh[iDir];

            std::vector<vec2cd> radPat(leaf->srcs.size(), vec2cd::Zero());

            int iSrc = 0;
            for (const auto& src : leaf->srcs)
                radPat[iSrc++] = toThPh * src->getRadAlongDir(center, kvec);

            leaf->radPats.push_back(radPat);
        }
    }
}

/* buildMpoleCoeffs()
 * (S2M) Build multipole coefficients from sources in this node  
 */
void FMM::Leaf::buildMpoleCoeffs() {
    if (isSrcless() || isRoot()) return;

    std::fill(coeffs.begin(), coeffs.end(), vec2cd::Zero());

    auto start = Clock::now();

    for (int iDir = 0; iDir < coeffs.size(); ++iDir) {
        vec2cd coeff = vec2cd::Zero();

        int iSrc = 0;
        for (const auto& src : srcs)
            coeff += (*lvec)[src->getIdx()] * radPats[iDir][iSrc++];

        coeffs[iDir] = coeff;
    }

    t.S2M += Clock::now() - start;
}

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void FMM::Leaf::buildLocalCoeffs() {
    if (isRoot()) return;

    auto start = Clock::now();
    buildMpoleToLocalCoeffs();
    t.M2L += Clock::now() - start;

    evalLeafIlistSols();

    start = Clock::now();
    if (!base->isRoot()) {
        localCoeffs = localCoeffs 
            + dynamic_cast<Stem*>(base)->getShiftedLocalCoeffs(branchIdx);
    }
    t.L2L += Clock::now() - start;
}

/* evalFarSols()
 * (L2T) Evaluate sols from local expansion due to far nodes
 */
//
void FMM::Leaf::evalFarSols() {
    if (isSrcless() || level <= 1) return;

    const auto [nth, nph] = angles[level].getNumAngles();

    int iObs = 0;
    for (const auto& obs : srcs) {
        size_t iDir = 0;
        cmplx intRad = 0;

        for (int ith = 0; ith < nth; ++ith) {
            for (int iph = 0; iph < nph; ++iph) {
                // Do the angular integration
                intRad += radPats[iDir][iObs].dot(localCoeffs[iDir]); // Hermitian dot!

                ++iDir;
            }
        }

        (*rvec)[obs->getIdx()] += Phys::C * wavenum * intRad;

        ++iObs;
    }
}
//

/*
void FMM::Leaf::evalFarSols() {
    if (isSrcless() || level <= 1) return;

    const auto [nth, nph] = angles[level].getNumAngles();

    const double phiWeight = 2.0*PI / static_cast<double>(nph); // TODO: static member

    size_t iObs = 0;
    for (const auto& obs : srcs) {
        size_t iDir = 0;
        cmplx intRad = 0;

        for (int ith = 0; ith < nth; ++ith) {
            const double weight = angles[level].weights[ith];

            for (int iph = 0; iph < nph; ++iph) {
                // Do the angular integration
                intRad += weight 
                    * radPats[iDir][iObs].dot(localCoeffs[iDir]); // Hermitian dot!

                ++iDir;
            }
        }

        (*rvec)[obs->getIdx()] += Phys::C * wavenum * phiWeight * intRad;

        ++iObs;
    }
}
*/

/* evalNearNonNborSols()
 * (M2T/S2T) Evaluate sols from mpole expansion due to list 3 nodes
 */
void FMM::Leaf::evalNearNonNborSols() {
    for (const auto& node : nearNonNbors)
        evalPairSols(node, nonNearRads[nonNearPairIdx++]);
    return;
}

/* (S2T) Evaluate sols at sources in this leaf due to sources in srcNode
 * and vice versa
 * srcNode : source node
 * rads    : precomputed radiation coefficients
 */
void FMM::Leaf::evalPairSols(const std::shared_ptr<Node> srcNode, const cmplxVec& rads) {

    const auto& srcSrcs = srcNode->getSrcs();

    const int nObs = srcs.size(), nSrcs = srcSrcs.size();

    cmplxVec solAtObss(nObs, 0.0);
    cmplxVec solAtSrcs(nSrcs, 0.0);

    int pairIdx = 0;
    for (size_t iObs = 0; iObs < nObs; ++iObs) {
        for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
            const auto obs = srcs[iObs], src = srcSrcs[iSrc];

            const cmplx rad = rads[pairIdx++];

            solAtObss[iObs] += (*lvec)[src->getIdx()] * rad;
            solAtSrcs[iSrc] += (*lvec)[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < nObs; ++n)
        (*rvec)[srcs[n]->getIdx()] += Phys::C * wavenum * solAtObss[n];

    for (int n = 0; n < nSrcs; ++n)
        (*rvec)[srcSrcs[n]->getIdx()] += Phys::C * wavenum * solAtSrcs[n];
}

/* evalSelfSols()
 * (S2T) Evaluate sols at sources in this node due to other sources
 * in this node
 */
void FMM::Leaf::evalSelfSols() {

    const int nSrcs = srcs.size();

    cmplxVec solAtObss(nSrcs, 0.0);

    // TODO: Handle self-interactions
    int pairIdx = 0;
    for (size_t iObs = 1; iObs < nSrcs; ++iObs) { // iObs = 0
        for (size_t iSrc = 0; iSrc < iObs; ++iSrc) { // iSrc <= iObs 
            auto obs = srcs[iObs], src = srcs[iSrc];

            const cmplx rad = selfRads[pairIdx++];

            solAtObss[iObs] += (*lvec)[src->getIdx()] * rad;
            solAtObss[iSrc] += (*lvec)[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < nSrcs; ++n)
        (*rvec)[srcs[n]->getIdx()] += Phys::C * wavenum * solAtObss[n];

}

/* evaluateSols()
 * Sum solutions at all sources in all leaves 
 */ 
void FMM::Leaf::evaluateSols() {
    auto start = Clock::now();
    for (const auto& leaf : leaves)
        leaf->evalFarSols();
    t.L2T += Clock::now() - start;

    for (const auto& [obsLeaf, srcLeaf] : nearPairs) {
        auto pairIdx = obsLeaf->leafPairIdx++;
        obsLeaf->evalPairSols(srcLeaf, obsLeaf->nearRads[pairIdx]);
    }

    for (const auto& leaf : leaves) {
        leaf->evalNearNonNborSols();

        leaf->evalSelfSols();

        leaf->leafPairIdx = 0;
        leaf->nonNearPairIdx = 0;
    }
}