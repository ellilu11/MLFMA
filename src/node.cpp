#include "node.h"

Config Node::config;
double Node::wavenum;
int Node::maxLevel;
std::vector<realVec> Node::phis;
std::vector<realVec> Node::thetas;
std::vector<realVec> Node::thetaWeights;
std::vector<int> Node::Ls;
Tables Node::tables;

void Node::setNodeParams(
    const Config& config_, const std::shared_ptr<Src>& Einc) {

    config = config_;
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
    nodeLeng(base == nullptr ? config.rootLeng : base->nodeLeng / 2.0),
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

    constexpr double EPS = 1.0E-9;

    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        const double nodeLeng = config.rootLeng / pow(2.0, lvl);

        // Use excess bandwidth formula
        const int L = ceil(
                (1.73*wavenum*nodeLeng +
                2.16*pow(config.digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));
            
        Ls.push_back(L);

        const int nph = 2*(L+1);
        realVec nodesPhi;
        for (int iph = 0; iph < nph; ++iph)
            nodesPhi.push_back(2.0*PI*iph/static_cast<double>(nph));
        phis.push_back(nodesPhi);

        const auto [nodes, weights] = Interp::gaussLegendre(L+1, EPS, 0.0, PI);
            // Interp::gaussLegendre(26, EPS, 0.0, PI);

        thetas.push_back(nodes);
        thetaWeights.push_back(weights);

        std::cout << "   (Lvl,Nth,Nph) = "
                  << "(" << lvl << "," << L+1 << "," << nph << ")\n";
    }

}

/* buildInteractionList()
 * Find interaction nodes
 */
void Node::buildInteractionList() {
    assert(!isRoot());
    assert(!nbors.empty());

    auto notContains = 
        [](const NodeVec& vec, const std::shared_ptr<Node>& val) {
        return std::find(vec.begin(), vec.end(), val) == vec.end();
    };

    for (const auto& baseNbor : base->nbors) {
        if (baseNbor->isSrcless()) continue; // TODO: double check

        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor);
            continue;
        }

        for (const auto& branch : baseNbor->branches)
            if (notContains(nbors, branch) && !branch->isSrcless()) // TODO: double check
                iList.push_back(branch);
    }

    assert(iList.size() <= pow(6, DIM) - pow(3, DIM));

    // std::cout << iList.size() << ' ';
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

    const auto order = config.interpOrder;

    const auto [nth, nph] = getNumAngles(level);
    localCoeffs.resize(nth*nph, vec3cd::Zero());

    const int nps = std::floor(q*(nth-1));

    for (const auto& node : iList) {
        const auto& mpoleCoeffs = node->getMpoleCoeffs();

        const auto& R = center - node->getCenter();
        const auto r = R.norm();
        const auto& rhat = R / r;
        
        const double normedDist = r / nodeLeng;

        size_t idx = 0;
        for (int ith = 0; ith < nth; ++ith) {

            for (int iph = 0; iph < nph; ++iph) {

                const auto& kvec = tables.kvec[level][idx];

                const double psi = kvec.dot(rhat) / wavenum;

                const int s = std::floor(nps * psi / PI);

                realVec psis;
                for (int ips = s+1-order; ips <= s+order; ++ips)
                    psis.push_back(PI*ips/static_cast<double>(nps));

                cmplx translCoeff = 0.0;

                for (int k = 0; k < 2*order; ++k) {
                    translCoeff +=
                        tables.transl[level].at(normedDist)[k]
                        * Interp::evalLagrangeBasis(psi,psis,k);
                        // * tables.interpPsi[level].at(psi)[k];
                }

                localCoeffs[idx] += translCoeff * mpoleCoeffs[idx];

                idx++;
            }
        }
    }
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

/* getFarSolsFromCoeffs(r)
 * Return sols at all sampled directions at distance r
 * from mpole coefficients of this node
 */
std::vector<vec3cd> Node::getFarSolsFromCoeffs(double r) {
    assert(!rwgs.empty());

    const cmplx C = -iu * c0 * wavenum * mu0
        * exp(iu*wavenum*r) / (4.0*PI*r);

    const auto [nth, nph] = getNumAngles(level);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto kvec = tables.kvec[level][idx];

            sols[idx] = 
                C * exp(-iu*kvec.dot(center)) 
                // * tables.matFromSph[level][idx]
                * coeffs[idx];

            idx++;
        }
    }

    return sols;
}

/* getFarSol()
 * Return sols at all sampled directions at distance r 
 * due to all RWGs in this node using farfield approximation
 */
std::vector<vec3cd> Node::getFarSols(double r) {

    assert(r >= 5.0 * config.rootLeng); // verify farfield condition

    const cmplx C = -iu * c0 * wavenum * mu0
        * exp(iu*wavenum*r) / (4.0*PI*r);

    const auto [nth, nph] = getNumAngles(level);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto& ImKK = tables.ImKK[level][idx];
            const auto& kvec = tables.kvec[level][idx];

            vec3cd dirCoeff = vec3cd::Zero();
            for (const auto& rwg : rwgs)
                dirCoeff += rwg->getRadAlongDir(center, kvec);

            sols[idx] = C * ImKK * dirCoeff;

            idx++;
        }
    }

    return sols;
}
