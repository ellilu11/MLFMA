#pragma once

#include "node.h"
#include "stem.h"

class Leaf;

using LeafVec = std::vector<std::shared_ptr<Leaf>>;
using LeafPair = std::pair<std::shared_ptr<Leaf>, std::shared_ptr<Leaf>>;

class Leaf final : public Node, public std::enable_shared_from_this<Leaf> {

public:
    Leaf(const SrcVec&, const int, Stem* const);

    std::shared_ptr<Node> getSelf() override { return shared_from_this(); }

    void pushToNearNonNbors(const std::shared_ptr<Node>& node) {
        nearNonNbors.push_back(node);
    }

    static int getNumLeaves() { return leaves.size(); }

    NodeVec getNearNonNbors() const { return nearNonNbors; }

    NodeVec getNearNbors() const { return nearNbors; }

    void evalFarSols();

    void evalNearNonNborSols();

    static std::vector<LeafPair> findNearNborPairs();

    static void evaluateSols();

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    void buildLocalCoeffs() override;

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << label << '\n';
    }

    static void printLeaves(std::ofstream& f) {
        for (const auto& leaf : leaves) {
            if (leaf->getSrcs().empty()) continue;

            f << leaf->getCenter() << " " << leaf->getLeng() << " " << leaf->getLevel() << '\n';
        }
    }

protected:
    static LeafVec leaves;

    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3
};
