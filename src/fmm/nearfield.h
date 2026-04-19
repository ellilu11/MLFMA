#pragma once

#include "node.h"

class FMM::Nearfield {

public:
    Nearfield(size_t);

    void evaluateSols();

    std::vector<NodePair> getNearPairs() const { return nearPairs; }

    std::vector<NodePair> getSelfPairs() const { return selfPairs; }

    void printNearMatrix(const std::string&) const;

private: 
    void findNodePairs();

    void buildTriPairs();

    size_t getNearCapacity();

    void buildNearMatrix();

public:
    sparseMat<cmplx> nearMat;

private:
    std::vector<NodePair> nearPairs;
    std::vector<NodePair> selfPairs;
};