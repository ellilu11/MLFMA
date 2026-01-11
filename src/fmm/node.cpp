#include "node.h"

Config FMM::Node::config;
double FMM::Node::wavenum;
std::vector<FMM::Angles> FMM::Node::angles;
std::vector<FMM::Tables> FMM::Node::tables;
std::vector<NodePair> FMM::Node::nonNearPairs;

std::shared_ptr<vecXcd> FMM::Node::lvec;
std::shared_ptr<vecXcd> FMM::Node::rvec;
std::shared_ptr<vecXcd> FMM::Node::currents;

void FMM::Node::initParams(
    const Config& config_,
    const std::shared_ptr<Excitation::PlaneWave>& Einc,
    int nsrcs)
{
    config = config_;
    wavenum = Einc->wavenum;

    lvec = std::make_shared<vecXcd>(vecXcd::Zero(nsrcs));
    rvec = std::make_shared<vecXcd>(vecXcd::Zero(nsrcs));
    currents = std::make_shared<vecXcd>(vecXcd::Zero(nsrcs)); // assume I = 0 initially

    Leaf::resetLeaves();
}

void FMM::Node::buildTables() { 
    std::cout << "   (Lvl,Nth,Nph) =\n";
    angles.reserve(maxLevel+1);
    for (int level = 0; level <= maxLevel; ++level)
        angles.emplace_back(level);

    Tables::buildDists();
    tables.reserve(maxLevel+1);
    for (int level = 0; level <= maxLevel; ++level)
        tables.emplace_back(level, maxLevel);
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
    const auto [nth, nph] = angles[level].getNumAngles();

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
        if (baseNbor->isSrcless()) continue;

        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor);
            continue;
        }

        for (const auto& branch : baseNbor->branches)
            if (notContains(nbors, branch) && !branch->isSrcless())
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

    const auto& transl = tables[level].transl;
    const size_t nDir = localCoeffs.size();

    for (const auto& node : iList) {
        const auto& dX = center - node->center;
        const auto& transl_dX = transl.at(dX/nodeLeng);

        for (int iDir = 0; iDir < nDir; ++iDir) // TODO: use Eigen::Array
            localCoeffs[iDir] += transl_dX[iDir] * node->coeffs[iDir];
    }

    // Apply integration weights
    const auto& angles_lvl = angles[level];
    const auto [nth, nph] = angles_lvl.getNumAngles();
    const double phiWeight = 2.0*PI / static_cast<double>(nph);

    size_t iDir = 0;
    for (int ith = 0; ith < nth; ++ith) 
        for (int iph = 0; iph < nph; ++iph)
            localCoeffs[iDir++] *= phiWeight * angles_lvl.weights[ith];
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
}

void FMM::Node::printFarSols(const std::string& fname) {
    namespace fs = std::filesystem;
    fs::path dir = "out/ff";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream farfile(dir/fname);
    farfile << std::setprecision(15) << std::scientific;

    const auto& angles_lvl = angles[level];
    const auto [nth, nph] = angles_lvl.getNumAngles();
    for (int iDir = 0; iDir < nth*nph; ++iDir) {
        const auto& krhat = angles_lvl.khat[iDir] * wavenum;

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += (*currents)[src->getIdx()] * src->getFarAlongDir(krhat);

        const vec3cd& far = Phys::C * wavenum * angles_lvl.ImRR[iDir] * dirFar;

        farfile << far << '\n';
    }

    // Also print out angles (coordinates of farsols)
    std::ofstream thfile(dir/"thetas.txt"), phfile(dir/"phis.txt");
    angles[level].printAngles(thfile, phfile);
}