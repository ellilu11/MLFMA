#pragma once

#include "node.h"

class Stem final : public Node, public std::enable_shared_from_this<Stem> {

public:
    Stem(const SrcVec&, const int, Stem* const);

    std::shared_ptr<Node> getSelf() override {
        return shared_from_this();
    }

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    const std::vector<vec3cd> getShiftedLocalCoeffs(int) const;

    const std::vector<vec3cd> getShiftedMpoleCoeffs(int) const; //

    static std::vector<vec3cd> getAnterpCoeffs(const std::vector<vec3cd>&, int);

    void buildLocalCoeffs() override;

    // ========== Test methods ==========
    void printNode(std::ofstream& f) override {
        f << center << " " << nodeLeng << " " << label << '\n';
        for (const auto& branch : branches)
            branch->printNode(f);
    }

    void printLocalCoeffs(std::ofstream& f) override {

        for (const auto& coeffs : localCoeffs)
            f << coeffs << '\n';

        /*
        for (const auto& branch : branches)
            branch->printLocalCoeffs(f);
            */
    }

    void testShiftedLocalCoeffs();

private:
    // TODO: Add branches
};
