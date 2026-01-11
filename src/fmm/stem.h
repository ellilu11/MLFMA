#pragma once

#include "node.h"

class FMM::Stem : public Node, public std::enable_shared_from_this<Stem> {

public:
    Stem(const SrcVec&, const int, Stem* const);

    void buildNeighbors() override;

    void buildLists() override;

    void resizeCoeffs() override;

    void buildMpoleCoeffs() override;

    std::vector<vec2cd> getShiftedLocalCoeffs(int) const;

    template <typename T>
    static void addInterpCoeffs(const std::vector<T>&, std::vector<T>&, int, int);

    template <typename T>
    static void addAnterpCoeffs(const std::vector<T>&, std::vector<T>&, int, int);

    void buildLocalCoeffs() override;

    std::shared_ptr<Node> getSelf() override { return shared_from_this(); }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << '\n';
    }

    static void tInterpPhi(int, int);

    static void tAnterpPhi(int, int);

    static void tInterpTheta(int, int);

    static void tAnterpTheta(int, int);

    template <typename T>
    static void tInterp(int, int);

    template <typename T>
    static void tAnterp(int, int);
};
