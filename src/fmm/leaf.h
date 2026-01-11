#pragma once

#include "node.h"
#include "stem.h"

class FMM::Leaf;

using LeafVec = std::vector<std::shared_ptr<FMM::Leaf>>;
using LeafPair = std::pair<std::shared_ptr<FMM::Leaf>, std::shared_ptr<FMM::Leaf>>;

class FMM::Leaf final : public Node, public std::enable_shared_from_this<Leaf> {

public:
    Leaf(const SrcVec&, const int, Stem* const);

    void buildLists() override;

    void resizeCoeffs() override;

    static void buildPairNearRads(const std::vector<NodePair>&, bool);

    static void buildNearRads();

    static void buildRadPats();

    void buildMpoleCoeffs() override;

    void buildLocalCoeffs() override;

    static void evaluateSols();

    std::shared_ptr<Node> getSelf() override { return shared_from_this(); }

    void pushToNearNonNbors(const std::shared_ptr<Node>& node) {
        nearNonNbors.push_back(node);
    }

    static int getNumLeaves() { return leaves.size(); }

    static void resetLeaves() { 
        leaves.clear(); 
        nearPairs.clear();
    }

    NodeVec getNearNonNbors() const { return nearNonNbors; }

    NodeVec getNearNbors() const { return nearNbors; }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << '\n';
    }

private:
    void buildNeighbors() override;

    static void findNearNborPairs();

    void evalFarSols();

    void evalNearNonNborSols();

    void evalPairSols(const std::shared_ptr<Node>, const cmplxVec&);

    void evalSelfSols();

    static LeafVec leaves;
    static std::vector<LeafPair> nearPairs;
    inline static size_t glSrcIdx = 0;

    std::vector<cmplxVec> nearRads;
    std::vector<cmplxVec> nonNearRads;
    cmplxVec selfRads;

    std::vector<std::vector<vec2cd>> radPats;
    std::array<std::vector<vec2cd>, 2> polarRadPats;

    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3

    size_t leafPairIdx; // index of near pairs with this leaf as observer
    size_t nonNearPairIdx; // index of non-near pairs with this leaf as observer
};
