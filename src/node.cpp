#include "node.h"

Config Node::config;
double Node::wavenum;
std::vector<realVec> Node::thetas;
std::vector<realVec> Node::thetaWeights;
std::vector<realVec> Node::phis;
std::vector<int> Node::Ls;

Tables Node::tables;
NodeVec Node::nodes;

void Node::setNodeParams(
    const Config& config_, const std::shared_ptr<PlaneWave>& Einc) 
{
    config = config_;
    wavenum = Einc->wavenum;
}

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
Node::Node(
    const SrcVec& srcs,
    const int branchIdx,
    Node* const base)
    : srcs(srcs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? config.rootLeng : base->nodeLeng / 2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? zeroVec :
        base->center + nodeLeng / 2.0 * Math::idx2pm(branchIdx))
{
    nodeIdx = numNodes++;
}

/* buildAngularSamples()
 * Compute theta and phi samples at each level
 */
void Node::buildAngularSamples() {

    constexpr double EPS_NR = 1.0E-9; // Newton-Raphson precision

    std::cout << "   (Lvl,Nth,Nph) =\n";

    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        const double nodeLeng = config.rootLeng / pow(2.0, lvl);

        // Use excess bandwidth formula
        const int tau = 
            ceil(
                (1.73*wavenum*nodeLeng +
                2.16*pow(config.digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));
                
        Ls.push_back(floor(0.50*tau)); // TODO: Find optimal formula

        // Construct thetas
        const int nth = tau+1;
        auto [nodes, weights] = Interp::gaussLegendre(nth, EPS_NR, 0.0, PI);

        // Absorb sin(theta) into weights
        std::transform(weights.begin(), weights.end(), nodes.begin(), weights.begin(),
            [](double weight, double theta) {
                return weight * sin(theta);
            });

        thetas.push_back(nodes);
        thetaWeights.push_back(weights);

        // Construct phis
        const int nph = 2*nth;
        realVec phis_lvl(nph);

        for (int iph = 0; iph < nph; ++iph)
            phis_lvl[iph] = 2.0*PI*iph/static_cast<double>(nph);

        phis.push_back(phis_lvl);

        std::cout << "   (" << lvl << "," << nth << "," << nph << ")\n";
    }
}

/* buildInteractionList()
 * Find interaction nodes
 */
void Node::buildInteractionList() {
    assert(!isRoot());
    assert(!nbors.empty());

    auto notContains = 
        [](const NodeVec& vec, const std::shared_ptr<Node> val) {
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

    const auto [nth, nph] = getNumAngles(level); 
    localCoeffs.resize(nth*nph, vec2cd::Zero()); // TODO: Allocate elsewhere

    if (iList.empty()) return;

    const int order = config.interpOrder;

    const int nps = std::floor(config.overInterp*(nth-1));

    for (const auto& node : iList) {

        const auto& mpoleCoeffs = node->coeffs;

        const auto& dX = center - node->center;
        const double r = dX.norm();
        const auto& rhat = dX / r;

        const auto& transl_dX = tables.transl[level].at(r/nodeLeng);

        for (int idx = 0; idx < nth*nph; ++idx) {
            const auto& khat = tables.khat[level][idx];

            const double psi = acos(khat.dot(rhat));

            cmplx translCoeff = 0.0;

            // psi LUT
            const auto [interpPsi, nearIdx] = tables.interpPsi[level].at(psi);

            for (int ips = nearIdx+1-order, k = 0; k < 2*order; ++ips, ++k) {

                const int ips_flipped = Math::flipIdxToRange(ips, nps); 

                translCoeff += transl_dX[ips_flipped] * interpPsi[k];
            }
            //

            localCoeffs[idx] += translCoeff * mpoleCoeffs[idx];
        }
    }
}

/* evalLeafIlistSols()
 * (S2L/S2T) Add contribution from list 4 to local coeffs
 */
void Node::evalLeafIlistSols() {
     for (const auto& node : leafIlist)
        evalPairSols(node);
     return;

    /* No psi LUT
    const int nearIdx = std::floor((nps-1) * psi / PI);

    realVec psis_;
    for (int ips = nearIdx+1-order; ips <= nearIdx+order; ++ips)
        psis_.push_back(PI*ips/static_cast<double>(nps-1));

    for (int ips = nearIdx+1-order, k = 0; k < 2*order; ++ips, ++k) {
        const int ips_flipped = Math::flipIdxToRange(ips, nps);

        translCoeff +=
            transls[ips_flipped]
            * Interp::evalLagrangeBasis(psi,psis_,k);
    }
    */
}

/* evalPairSols(srcNode)
 * (S2T) Evaluate sols at sources in this node due to sources in srcNode
 * and vice versa
 * srcNode : source node
 */
void Node::evalPairSols(const std::shared_ptr<Node> srcNode) {

    const int numObss = srcs.size(), numSrcs = srcNode->srcs.size();

    cmplxVec solAtObss(numObss, 0.0);
    cmplxVec solAtSrcs(numSrcs, 0.0);

    for (size_t obsIdx = 0; obsIdx < numObss; ++obsIdx) {
        for (size_t srcIdx = 0; srcIdx < numSrcs; ++srcIdx) {
            const auto obs = srcs[obsIdx], src = srcNode->srcs[srcIdx];

            const cmplx rad = obs->getIntegratedRad(src);

            solAtObss[obsIdx] += src->getCurrent() * rad;
            solAtSrcs[srcIdx] += obs->getCurrent() * rad;
        }
    }

    for (int n = 0; n < numObss; ++n)
        srcs[n]->addToSol(Phys::C * wavenum * solAtObss[n]);

    for (int n = 0; n < numSrcs; ++n)
        (srcNode->srcs[n])->addToSol(Phys::C * wavenum * solAtSrcs[n]);
}

/* evalSelfSols()
 * (S2T) Evaluate sols at sources in this node due to other sources
 * in this node
 */
void Node::evalSelfSols() {

    const int numSrcs = srcs.size();

    cmplxVec sols(numSrcs, 0.0);

    for (size_t obsIdx = 1; obsIdx < numSrcs; ++obsIdx) { // obsIdx = 0 to include self-interactions
        for (size_t srcIdx = 0; srcIdx < obsIdx; ++srcIdx) {
            auto obs = srcs[obsIdx], src = srcs[srcIdx];

            const cmplx rad = obs->getIntegratedRad(src);

            sols[obsIdx] += src->getCurrent() * rad;
            sols[srcIdx] += obs->getCurrent() * rad;
        }
    }

    for (int n = 0; n < numSrcs; ++n)
        srcs[n]->addToSol(Phys::C * wavenum * sols[n]);
}


/*void Node::evalPairSols(const std::shared_ptr<Node> srcNode) {

    assert(getSelf() != srcNode);

    const auto& srcSrcs = srcNode->srcs;

    for (const auto& obs : srcs) {
        cmplx sol = 0;

        for (const auto& src : srcSrcs)
            sol += src->getCurrent() * obs->getIntegratedRad(src);

        obs->addToSol(C * wavenum * sol);
    }
}*/

void Node::evalSelfSolsSlow() {

    for (const auto& obs : srcs) {
        cmplx sol = 0;

        for (const auto& src : srcs) {
            if (src == obs) continue; // TODO: Use analytic expression

            sol += src->getCurrent() * obs->getIntegratedRad(src);

        }

        obs->addToSol(Phys::C * wavenum * sol);
    }
}

/* getFarSol()
 * Return sols at all sampled directions at distance r 
 * due to all sources in this node using farfield approximation
 */
/*std::vector<vec3cd> Node::getFarSols(double r) {

    assert(r >= 5.0 * config.rootLeng); // verify farfield condition

    const cmplx kernel = C * wavenum * exp(iu*wavenum*r) / r;

    const auto [nth, nph] = getNumAngles(level);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto& ImKK = tables.ImKK[level][idx];
            const auto& kvec = tables.khat[level][idx] * wavenum;

            vec3cd dirCoeff = vec3cd::Zero();
            for (const auto& src : srcs)
                dirCoeff += src->getRadAlongDir(center, kvec);

            sols[idx] = kernel * ImKK * dirCoeff;

            idx++;
        }
    }

    return sols;
}*/