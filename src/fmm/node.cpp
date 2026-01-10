#include "node.h"

Config FMM::Node::config;
double FMM::Node::wavenum;
FMM::Angles FMM::Node::angles;
FMM::Tables FMM::Node::tables;
std::vector<NodePair> FMM::Node::nonNearPairs;

std::shared_ptr<vecXcd> FMM::Node::lvec;
std::shared_ptr<vecXcd> FMM::Node::rvec;
std::shared_ptr<vecXcd> FMM::Node::currents;

void FMM::Node::initParams(
    const Config& config_,
    const std::shared_ptr<Excitation::PlaneWave>& Einc)
{
    config = config_;
    wavenum = Einc->wavenum;
}

void FMM::Node::buildTables() { 
    angles = Angles(wavenum, config.rootLeng, config.digits, maxLevel);
    tables = Tables(angles, config, wavenum, maxLevel); 
}

void FMM::Node::linkStates(const std::unique_ptr<Solver>& solver) {
    lvec = std::move(solver->getLvec());
    rvec = std::move(solver->getRvec());
    currents = std::move(solver->getCurrents());
}

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
FMM::Node::Node(
    const SrcVec& srcs,
    const int branchIdx,
    Node* const base)
    : srcs(srcs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? config.rootLeng : base->nodeLeng/2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? zeroVec :
        base->center + nodeLeng/2.0 * Math::idx2pm(branchIdx))
{
    ++numNodes;
}

void FMM::Node::resizeCoeffs() {
    const auto [nth, nph] = angles.getNumAngles(level);

    coeffs.resize(nth*nph, vec2cd::Zero());
    localCoeffs.resize(nth*nph, vec2cd::Zero());
}

/* buildInteractionList()
 * Find interaction nodes
 */
void FMM::Node::buildInteractionList() {
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
}

/* pushSelfToNearNonNbors()
 * Add this node to list 3 of leaf.
 * (if leaf is in list 4 of self, self is in list 3 of leaf) 
 */
void FMM::Node::pushSelfToNearNonNbors() {
    if (leafIlist.empty()) return;

    for (const auto& node : leafIlist) {
        auto leaf = dynamic_pointer_cast<Leaf>(node);

        leaf->pushToNearNonNbors(getSelf());
        nonNearPairs.emplace_back(leaf, getSelf()); // record list4-list3 pair
    }
}

/* buildMpoleToLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 */
void FMM::Node::buildMpoleToLocalCoeffs() {
    if (iList.empty()) return;

    std::fill(localCoeffs.begin(), localCoeffs.end(), vec2cd::Zero());

    const size_t numAngles = localCoeffs.size();

    for (const auto& node : iList) {
        const auto& dX = center - node->center;
        const auto& transl_dX = tables.transl[level].at(dX/nodeLeng);

        for (int idx = 0; idx < numAngles; ++idx) // TODO: use Eigen::Array
            localCoeffs[idx] += transl_dX[idx] * node->coeffs[idx];
    }

    // Apply integration weights
    const auto [nth, nph] = angles.getNumAngles(level);
    const double phiWeight = 2.0*PI / static_cast<double>(nph);
    size_t dirIdx = 0;
    for (int ith = 0; ith < nth; ++ith) {
        const double weight = phiWeight * angles.thetaWeights[level][ith];

        for (int iph = 0; iph < nph; ++iph)
            localCoeffs[dirIdx++] *= weight;
    }
    //
}

/* evalLeafIlistSols()
 * (S2L/S2T) Add contribution from list 4 to local coeffs
 */
void FMM::Node::evalLeafIlistSols() {
    // Do nothing! Contribution from list 4 node is 
    // to be evaluated by Leaf::evalNearNonNborSols()
    //for (const auto& node : leafIlist)
    //    evalPairSols(node, nonNearRads[nonNearPairIdx++]);
    //return;

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

void FMM::Node::printFarSols(const std::string& fname) {
    // assert(r >= 5.0 * config.rootLeng); // verify farfield condition
    // const cmplx kernel = Phys::C * wavenum * exp(iu*wavenum*r) / r;

    namespace fs = std::filesystem;
    fs::path dir = "out/ff";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream farfile(dir/fname);
    farfile << std::setprecision(15) << std::scientific;

    const auto [nth, nph] = angles.getNumAngles(level);
    for (int dirIdx = 0; dirIdx < nth*nph; ++dirIdx) {
        const auto& krhat = tables.khat[level][dirIdx] * wavenum;

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += (*currents)[src->getIdx()] * src->getFarAlongDir(krhat);

        const vec3cd& far = Phys::C * wavenum * tables.ImRR[level][dirIdx] * dirFar;

        farfile << far << '\n';
    }
}

void FMM::Node::printAngles() {
    std::filesystem::path dir = "out/ff";
    std::ofstream thfile(dir/"thetas.txt"), phfile(dir/"phis.txt");

    for (const auto& theta : angles.thetas[level])
        thfile << theta << '\n';

    for (const auto& phi : angles.phis[level])
        phfile << phi << '\n';
}