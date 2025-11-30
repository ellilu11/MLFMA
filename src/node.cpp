#include "node.h"

int Node::order;
int Node::orderExp;
int Node::maxNodeSrcs;
int Node::maxLevel;
double Node::rootLeng;
double Node::k;
int Node::L;
std::vector<realVec> Node::thetas;
std::vector<realVec> Node::thetaWeights;
Tables Node::tables;

void Node::setNodeParams(
    const Config& config, const std::shared_ptr<Src>& Einc) {

    order = config.order; // ceil(-std::log(config.EPS) / std::log(2));
    orderExp = [&]() -> std::size_t {
        switch (config.prec) {
            case Precision::LOW:    return 8;
            case Precision::MEDIUM: return 17;
            case Precision::HIGH:   return 26;
        }
        }();
    maxNodeSrcs = config.maxNodeSrcs;
    rootLeng = config.rootLeng; // TODO: define from max l_infty norm of all rwg centers
    k = Einc->k;
}

/* setThetaSamples(config)
 * Compute L and theta samples at each level
 */
void Node::setThetaSamples() {

    for (int lvl = maxLevel; lvl >= 0; --lvl) {
        const double nodeLeng = rootLeng / pow(2.0, lvl);

        // Use excess bandwidth formula
        L = ceil(1.73*k*nodeLeng +
            2.16*pow(orderExp, 2.0/3.0)*pow(k*nodeLeng, 1.0/3.0));

        const auto [nodes, weights] = Math::gaussLegendreTheta(L+1, 1.0E-9);

        thetas[lvl] = nodes;
        thetaWeights[lvl] = weights;
    }
}

void Node::buildTables(const Config& config) {
    tables = Tables(maxLevel,k,thetas);
    // assert(orderExp == tables.quadCoeffs_.size());
}

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
Node::Node(
    const RWGVec& rwgs,
    const int branchIdx,
    Node* const base)
    : rwgs(rwgs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? rootLeng : base->nodeLeng / 2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? zeroVec :
        base->center + nodeLeng / 2.0 * Math::idx2pm(branchIdx)),
    label(0)
{
    // for (int l = 0; l <= order; ++l) 
    //    localCoeffs.push_back(vec2cd::Zero(2*l+1));

    maxLevel = level;
    numNodes++;
}

/* buildInteractionList()
 * Find interaction nodes, and assign each to a dirlist.
 */
void Node::buildInteractionList() {
    assert(!isRoot());
    assert(!nbors.empty());

    auto notContains = [](const NodeVec& vec, const std::shared_ptr<Node>& val) {
        return std::find(vec.begin(), vec.end(), val) == vec.end();
    };

    NodeVec iList;
    for (const auto& baseNbor : base->nbors) {
        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor);
            continue;
        }
        for (const auto& branch : baseNbor->branches)
            if (notContains(nbors, branch))
                iList.push_back(branch);
    }

    assert(iList.size() <= pow(6, DIM) - pow(3, DIM));

    // pick minDist \in (nodeLeng, 2.0*nodeLeng) to avoid rounding errors
    const double minDist = 1.5 * nodeLeng;

}

/* pushSelfToNearNonNbors()
 * Add this node to list 3 of leaf.
 * (if leaf is in list 4 of self, self is in list 3 of leaf) 
 */
void Node::pushSelfToNearNonNbors() {
    if (leafIlist.empty()) return;

    for (const auto& node : leafIlist) {
        auto leaf = dynamic_pointer_cast<Leaf>(node);
        leaf->pushToNearNonNbors(getSelf()); // call shared_from_this()
    }
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
/*std::vector<vecXcd> Node::getShiftedLocalCoeffs(const int branchIdx) const {

}*/

/* evalLeafIlistSols()
 * (P2L) Add contribution from list 4 to local coeffs
 */
void Node::evalLeafIlistSols() {
  
}

/* evalPairSols(srcNode)
 * (P2P) Evaluate sols at particles in this node due to particles in srcNode
 * and vice versa 
 * srcNode : source node
 */
void Node::evalPairSols(const std::shared_ptr<Node>& srcNode) {
    const int numObss = rwgs.size(), numSrcs = srcNode->rwgs.size();
}

/* evalSelfSols()
 * (P2P) Evaluate sols at all particles in this node due to all other particles
 * in this node
 */
void Node::evalSelfSols() {
    const int numRWG = rwgs.size();

}

