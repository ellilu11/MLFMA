#pragma once

#include "node.h"

class Stem : public Node, public std::enable_shared_from_this<Stem> {

public:
    Stem(const SrcVec&, const int, Stem* const);

    std::shared_ptr<Node> getSelf() override { return shared_from_this(); }

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    std::vector<vec2cd> getShiftedLocalCoeffs(int) const;

    template <typename T>
    static void addInterpCoeffs(const std::vector<T>&, std::vector<T>&, int, int);

    template <typename T>
    static void addAnterpCoeffs(const std::vector<T>&, std::vector<T>&, int, int);

    void buildLocalCoeffs() override;

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << '\n';
    }

    static void tInterpPhi(int);

    static void tAnterpPhi(int);

    static void tInterp(int);

    static void tAnterp(int);
};
