#pragma once

#include "node.h"
#include "stem.h"

class Leaf;

using LeafVec = std::vector<std::shared_ptr<Leaf>>;
using LeafPair = std::pair<std::shared_ptr<Leaf>, std::shared_ptr<Leaf>>;

class Leaf final : public Node, public std::enable_shared_from_this<Leaf> {

public:
    Leaf(
        const RWGVec&,
        const int,
        Stem* const);

    std::shared_ptr<Node> getSelf() override {
        return shared_from_this();
    }

    void pushToNearNonNbors(const std::shared_ptr<Node>& node) {
        nearNonNbors.push_back(node);
    }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << label << '\n';
    }

    static const int getNumLeaves() { return leaves.size(); }

    const NodeVec getNearNonNbors() const { return nearNonNbors; }

    const NodeVec getNearNbors() const { return nearNbors; }

    void evalFarSols();

    void evalNearNonNborSols();

    static std::vector<LeafPair> findNearNborPairs();

    static void evaluateSols();

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    void buildLocalCoeffs() override;

    // ========== Test methods ==========
    vec2cd getLeafSols(const vec3d);

    static void testFarFieldFromLeaves(const std::vector<vec3d>&);

private:
    static LeafVec leaves;

    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3
};
