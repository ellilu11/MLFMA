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
        base->center + nodeLeng / 2.0 * Math::idx2pm(branchIdx)),
    label(0)
{
    nodeIdx = numNodes++;
}

/* buildAngularSamples()
 * Compute theta and phi samples at each level
 */
void Node::buildAngularSamples() {

    constexpr double EPS_NR = 1.0E-9; // Newton-Raphson precision

    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        const double nodeLeng = config.rootLeng / pow(2.0, lvl);

        // Use excess bandwidth formula
        const int tau = // 15; 
            ceil(
                (1.73*wavenum*nodeLeng +
                2.16*pow(config.digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));
            
        // TODO: Find an optimal formula, possibly explicitly depending on level
        Ls.push_back(tau/2 - 1); 

        // Construct thetas
        const int nth = tau+1;
        const auto [nodes, weights] = Interp::gaussLegendre(nth, EPS_NR, 0.0, PI);

        thetas.push_back(nodes);
        thetaWeights.push_back(weights);

        // Construct phis
        const int nph = 2*nth;
        realVec phis_lvl;

        for (int iph = 0; iph < nph; ++iph)
            phis_lvl.push_back(2.0*PI*iph/static_cast<double>(nph));

        phis.push_back(phis_lvl);

        std::cout << "   (Lvl,Nth,Nph) = "
                  << "(" << lvl << "," << nth << "," << nph << ")\n";
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
    localCoeffs.resize(nth*nph, vec3cd::Zero());

    if (iList.empty()) return;

    const int order = config.interpOrder;

    const int nps = std::floor(Q*(nth-1));

    for (const auto& node : iList) {

        const auto& mpoleCoeffs = node->coeffs;

        const auto& dR = center - node->center;
        const double r = dR.norm();
        const auto& rhat = dR / r;

        const auto& translVec = tables.transl[level].at(r / nodeLeng);

        size_t idx = 0;
        for (int ith = 0; ith < nth; ++ith) {

            for (int iph = 0; iph < nph; ++iph) {

                const auto& khat = tables.khat[level][idx];

                const double psi = acos(khat.dot(rhat));

                cmplx translCoeff = 0.0;

                /* No psi LUT
                const int s = std::floor((nps-1) * psi / PI);

                realVec psis_;
                for (int ips = s+1-order; ips <= s+order; ++ips)
                    psis_.push_back(PI*ips/static_cast<double>(nps-1));

                for (int ips = s+1-order, k = 0; k < 2*order; ++ips, ++k) {
                    const int ips_flipped = Math::flipIdxToRange(ips, nps);

                    translCoeff +=
                        translVec[ips_flipped]
                        * Interp::evalLagrangeBasis(psi,psis_,k);
                }
                */

                // psi LUT
                const auto& interpVec = tables.interpPsi[level].at(psi);
                const int s = tables.ssps[level].at(psi);

                for (int ips = s+1-order, k = 0; k < 2*order; ++ips, ++k) {

                    const int ips_flipped = Math::flipIdxToRange(ips, nps); 

                    translCoeff +=
                        translVec[ips_flipped]
                        * interpVec[k];
                }
                //

                localCoeffs[idx] += translCoeff * mpoleCoeffs[idx];

                idx++;
            }
        }
    }
}


// No psi interpolation
/*void Node::buildMpoleToLocalCoeffs() {
    if (iList.empty()) return;

    const int order = config.interpOrder;

    const auto [nth, nph] = getNumAngles(level);
    localCoeffs.resize(nth*nph, vec3cd::Zero());

    const int L = Ls[level];

    // for (const auto& node : iList) {
        const auto& node = iList[0];

        const auto& mpoleCoeffs = node->coeffs;

        const auto& dR = center - node->center;
        const double r = dR.norm();
        const auto& rhat = dR / r;

        const double kr = wavenum*r;
        vecXcd radialCoeffs(L+1);

        for (int l = 0; l <= L; ++l) {
            radialCoeffs[l] =
                iu * wavenum / (4.0*PI)
                * pow(iu,l) * (2.0*l+1.0)
                * Math::sphericalHankel1(kr, l);
        }

        size_t idx = 0;
        for (int ith = 0; ith < nth; ++ith) {

            for (int iph = 0; iph < nph; ++iph) {

                const auto& khat = tables.khat[level][idx];

                const double xi = khat.dot(rhat);

                cmplx translCoeff = 0.0;

                for (int l = 0; l <= L; ++l) {
                    translCoeff +=
                        radialCoeffs[l]
                        * Math::legendreP(xi, l).first
                        ;
                }

                localCoeffs[idx] += translCoeff * mpoleCoeffs[idx];

                idx++;
            }
        }
    // } // for (const auto& node : iList)
}*/

/* evalLeafIlistSols()
 * (S2L) Add contribution from list 4 to local coeffs
 */
void Node::evalLeafIlistSols() {
     for (const auto& node : leafIlist)
        evalPairSols(node);
     return;
}

/* evalPairSols(srcNode)
 * (S2T) Evaluate sols at sources in this node due to sources in srcNode
 * srcNode : source node
 */
void Node::evalPairSols(const std::shared_ptr<Node> srcNode) {

    assert(getSelf() != srcNode);

    const auto& srcSrcs = srcNode->srcs;

    for (const auto& obs : srcs) {
        cmplx sol = 0;

        for (const auto& src : srcSrcs)
             sol += obs->getIntegratedRad(src);

        obs->addToSol(C * wavenum * sol);
    }
}

/* evalSelfSols()
 * (S2T) Evaluate sols at sources in this node due to other sources
 * in this node
 */
void Node::evalSelfSols() {

    for (const auto& obs : srcs) {
        cmplx sol = 0;

        for (const auto& src : srcs) {
            if (src == obs) continue; // TODO: Use analytic expression

            sol += obs->getIntegratedRad(src);

        }

        obs->addToSol(C * wavenum * sol);
    }
}

/* getFarSol()
 * Return sols at all sampled directions at distance r 
 * due to all sources in this node using farfield approximation
 */
std::vector<vec3cd> Node::getFarSols(double r) {

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
}