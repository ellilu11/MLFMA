#pragma once

#include "node.h"

class Stem final : public Node, public std::enable_shared_from_this<Stem> {

public:
    Stem(
        const RWGVec&,
        const int,
        Stem* const);

    std::shared_ptr<Node> getSelf() override {
        return shared_from_this();
    }

    void printNode(std::ofstream& f) override {
        f << center << " " << nodeLeng << " " << label << '\n';
        for (const auto& branch : branches)
            branch->printNode(f);
    }

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    void propagateExpCoeffs() override;

    void buildLocalCoeffs() override;

};
