#include "node.h"

int Node::order;
int Node::prec;
int Node::maxNodeSrcs;
int Node::maxLevel;
double Node::rootLeng;
double Node::wavenum;
std::vector<realVec> Node::phis;
std::vector<realVec> Node::thetas;
std::vector<realVec> Node::thetaWeights;
std::vector<int> Node::Ls;
Tables Node::tables;

void Node::setNodeParams(
    const Config& config, const std::shared_ptr<Src>& Einc) {

    order = config.order; // ceil(-std::log(config.EPS) / std::log(2));
    prec = [&]() -> std::size_t {
        switch (config.prec) {
            case Precision::LOW:    return 3;
            case Precision::MEDIUM: return 6;
            case Precision::HIGH:   return 9;
        }
        }();
    maxNodeSrcs = config.maxNodeSrcs;
    rootLeng = config.rootLeng;
    wavenum = Einc->wavenum;
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
    numNodes++;
}

/* buildAngularSamples()
 * Compute theta and phi samples at each level
 */
void Node::buildAngularSamples() {

    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        const double nodeLeng = rootLeng / pow(2.0, lvl);

        // Use excess bandwidth formula
        const int L = ceil(1.73*wavenum*nodeLeng +
            2.16*pow(prec, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0));

        Ls.push_back(L);

        const int nph = 2*(L+1);
        realVec nodesPhi;
        for (int iph = 0; iph < nph; ++iph)
            nodesPhi.push_back(2.0*PI*iph/static_cast<double>(nph));
        phis.push_back(nodesPhi);

        const auto [nodes, weights] = Interp::gaussLegendre(L+1, 1.0E-9, 0.0, PI);

        thetas.push_back(nodes);
        thetaWeights.push_back(weights);
    }
}

void Node::buildTables(const Config& config) {
    tables = Tables(
        maxLevel, config.order, thetas, phis, Ls, wavenum, rootLeng);
}

/* buildInteractionList()
 * Find interaction nodes
 */
void Node::buildInteractionList() {
    assert(!isRoot());
    assert(!nbors.empty());

    auto notContains = [](const NodeVec& vec, const std::shared_ptr<Node>& val) {
        return std::find(vec.begin(), vec.end(), val) == vec.end();
    };

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

/* buildMpoleToLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 */
void Node::buildMpoleToLocalCoeffs() {
    if (isRoot()) return;

    const auto [nth, nph] = getNumAngles(level);
    localCoeffs.resize(nth*nph, vec2cd::Zero());

    /*for (const auto& node : iList) {
        const auto mpoleCoeffs = node->getMpoleCoeffs();

        const auto R = center - node->getCenter();
        const auto r = R.norm();
        const auto rhat = R / r;
        
        const auto iDist = Math::getIdxOfVal(R, tables.normedDists);

        const int t = 0; // std::floor()

        size_t idx = 0;
        for (int ith = 0; ith < nth; ++ith) {
            for (int iph = 0; iph < nph; ++iph) {

                const vec3d kvec = tables.kvec[level][idx];

                const double psi = kvec.dot(rhat) / kvec.norm();

                auto lagrangeCoeff = Interp::evalLagrangeBasis(psi, psis, i);

                double translCoeff = 0.0;

                for (int i = t+1-order; i <= t+order; ++i) {
                    translCoeff +=
                        tables.transl[level][iDist][i] *
                        lagrangeCoeff;
                }

                localCoeffs[idx] += translCoeff * mpoleCoeffs[idx];

                idx++;
            }
        }
    }*/
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
std::vector<vec2cd> Node::getShiftedLocalCoeffs(const int branchIdx) const {
    std::vector<vec2cd> shiftedLocalCoeffs;

    return shiftedLocalCoeffs;
}

/* evalLeafIlistSols()
 * (S2L) Add contribution from list 4 to local coeffs
 */
void Node::evalLeafIlistSols() {
  
}

/* evalPairSols(srcNode)
 * (S2T) Evaluate sols at RWGs in this node due to RWGs in srcNode
 * and vice versa 
 * srcNode : source node
 */
void Node::evalPairSols(const std::shared_ptr<Node>& srcNode) {
    const int numObss = rwgs.size(), numSrcs = srcNode->rwgs.size();
}

/* evalSelfSols()
 * (S2T) Evaluate sols at all RWGs in this node due to all other RWGs
 * in this node
 */
void Node::evalSelfSols() {
    const int numRWGs = rwgs.size();

}

/* getFarSols()
 * Return sols at spherical point R due to all RWGs in this node
 * using farfield approximation
 */
vec2cd Node::getFarSols(const vec3d R) {
    auto r = R[0];
    auto rhat = Math::fromSph(R) / r;
    auto ImRR = Math::IminusRR(R[1], R[2]);

    vec3cd fvec = vec3cd::Zero();

    for (const auto& rwg : rwgs) {
        vec3cd rwgCoeff = vec3cd::Zero();

        auto triPlus = rwg->getTriPlus();
        auto [nodesPlus, weightPlus] = triPlus->getQuads();
        for (const auto& quadNode : nodesPlus)
            rwgCoeff += weightPlus * Math::expI(-wavenum*rhat.dot(quadNode)) 
                        * (rwg->getVplus() - quadNode);

        auto triMinus = rwg->getTriMinus();
        auto [nodesMinus, weightMinus] = triMinus->getQuads();
        for (const auto& quadNode : nodesMinus)
            rwgCoeff += weightMinus * Math::expI(-wavenum*rhat.dot(quadNode)) 
                        * (quadNode - rwg->getVminus());
        
        fvec += rwg->getCurrent() * rwg->getLeng() * rwgCoeff;
    }

    return -iu * wavenum 
            * Math::expI(wavenum*r) / (4.0*PI*r) 
            * Math::matToThPh(R[1], R[2]) * (ImRR * fvec);
}
