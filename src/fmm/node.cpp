#include "node.h"
#include "../mesh/rwg.h"

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
FMM::Node::Node(
    const SrcVec& srcs,
    const int branchIdx,
    Node* const base,
    bool buildLeaf)
    : srcs(srcs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? Mesh::rootLeng : base->nodeLeng/2.0),
    center(base == nullptr ? Mesh::rootCenter :
        base->center + nodeLeng/2.0 * Math::idx2pm(branchIdx)),
    lvl(base == nullptr ? 0 : base->lvl + 1)
{
    if (isRoot()) std::cout << " Building FMM tree...\n";

    if (buildLeaf) {
        // Assign contiguous global indices to srcs in this leaf node
        // TODO: Use Morton ordering to improve spatial locality of srcs in leaf nodes
        for (const auto& src : srcs) src->setIdx(glSrcIdx++);

        // Update maxLevel if needed
        maxLevel = std::max(lvl, maxLevel);
    } else subdivide();

    // Find unique triangles of RWGs in this node
    // TODO: Only call for leaves and non-near stems
    findTris();
}

/* subdivide()
 * Subdivide node into 8 branches and assign srcs to branches
 */
void FMM::Node::subdivide() {
    auto bools2Idx = [](const vec3i& bools) {
        return bools[0] + 2*bools[1] + 4*bools[2];
    };

    // Assign every src to a branch based on src center relative to node center
    std::array<SrcVec, 8> branchSrcs;
    for (const auto& src : srcs) {
        size_t idx = bools2Idx(src->getCenter() > center);
        branchSrcs[idx].push_back(src);
    }

    // Construct branch nodes
    branches.resize(8);
    for (size_t k = 0; k < 8; ++k) {
        bool buildLeaf = branchSrcs[k].size() <= config.maxNodeSrcs;
        branches[k] = std::make_shared<Node>(branchSrcs[k], k, this, buildLeaf);
    }
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 * If node is leaf, find all neighbor leaves of equal or lesser size (list 1)
 */
void FMM::Node::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (!nbor) continue; // continue if nbor is nullptr
        
        nbors.push_back(nbor);

        if (!isLeaf()) continue;
        auto nbors = getNeighborsLeqSize(nbor, dir);
        nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
    }

    assert(nbors.size() <= numDir);
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

        if (baseNbor->isLeaf() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor); // TODO: Subdivide baseNbor instead of adding to leafIlist
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
        node->pushToNearNonNbors(shared_from_this());
        glNonNearPairs.emplace_back(node, shared_from_this()); // record list4-list3 pair
    }
}

/* postProcess()
 * Find neighbor and interaction lists.
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void FMM::Node::postProcess() {
    // If node has sources, add to global list of nodes/leaves
    if (!isSrcless()) {
        glNodes.push_back(shared_from_this());
        glNodesByLvl[lvl].push_back(shared_from_this());
        if (isLeaf()) glLeaves.push_back(shared_from_this());
    }
 
    // Find neighbor and interaction lists, and add self as near non-neighbor
    if (!isRoot()) {
        buildNeighbors();
        buildInteractionList();
        pushSelfToNearNonNbors();
    }

    // Recurse into branches
    for (const auto& branch : branches)
        branch->postProcess();

    if (isRoot()) {
        std::cout << "   # Nodes:         " << glNodes.size() << '\n';
        std::cout << "   # Leaves:        " << glLeaves.size() << '\n';
        std::cout << "   Max node level:  " << maxLevel << "\n\n";
    }
}

/* findTris()
 * Find all unique triangles of RWGs in this node
 */
void FMM::Node::findTris() {
    if (isSrcless()) return;

    for (const auto& src : srcs) {
        if (!src->isSrcType<Mesh::RWG>()) continue;
        const auto rwg = dynamic_pointer_cast<Mesh::RWG>(src);

        auto [iTri0, iTri1] = rwg->getTrisIdx();
        iTris.push_back(iTri0);
        iTris.push_back(iTri1);
    }

    std::sort(iTris.begin(), iTris.end());
    iTris.erase(std::unique(iTris.begin(), iTris.end()), iTris.end());
}
