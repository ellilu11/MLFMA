#pragma once

#include "node.h"

class Stem : public Node, public std::enable_shared_from_this<Stem> {

public:
    Stem(const SrcVec&, const int, Stem* const);

    std::shared_ptr<Node> getSelf() override { return shared_from_this(); }

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    std::vector<vec3cd> getShiftedLocalCoeffs(int) const;

    static void addInterpCoeffs(const std::vector<vec3cd>&, std::vector<vec3cd>&, int, int);

    // static void addAnterpCoeffs(const std::vector<vec3cd>&, std::vector<vec3cd>&, int);

    void buildLocalCoeffs() override;

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << label << '\n';
    }

    static void testInvInterp(int);

// protected:
    // TODO: Add branches
};
