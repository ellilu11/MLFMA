#pragma once

#include "node.h"

class Stem : public Node, public std::enable_shared_from_this<Stem> {

public:
    Stem(const SrcVec&, const int, Stem* const);

    std::shared_ptr<Node> getSelf() override { return shared_from_this(); }

    void buildNeighbors() override;

    void initNode() override;

    void buildMpoleCoeffs() override;

    std::vector<vec2cd> getShiftedLocalCoeffs(int) const;

    static void addInterpCoeffs(const std::vector<vec2cd>&, std::vector<vec2cd>&, int, int);

    // static void addAnterpCoeffs(const std::vector<vec2cd>&, std::vector<vec2cd>&, int);

    void buildLocalCoeffs() override;

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << '\n';
    }

    static void testInvInterp(int);
};
