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

    void buildNeighbors() override;

    void buildLists() override;

    static void findNearNborPairs();

    static void buildNearRads();

    static void buildRadPats();

    void buildMpoleCoeffs() override;

    void buildLocalCoeffs() override;

    void evalFarSols();

    void evalNearNonNborSols();

    void evalPairSols(const std::shared_ptr<Node>);

    void evalSelfSols();

    static void evaluateSols();

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << '\n';
    }

protected:
    static LeafVec leaves;
    static std::vector<LeafPair> nearPairs;

    cmplxVec nearRads;
    cmplxVec selfRads;

    std::vector<std::vector<vec2cd>> radPats;

    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3
};
