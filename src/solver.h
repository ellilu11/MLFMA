#pragma once

#include <filesystem>
#include "types.h"
#include "sources/source.h"

class Node;

class Solver {

public :
    Solver(SrcVec& srcs, std::shared_ptr<Node>, int, double);

    void updateSols(int);

    void iterateArnoldi(int);

    void applyGivensRotation(vecXcd&, vecXcd&, vecXcd&, int);

    void updateGvec(vecXcd&, vecXcd&, int);

    void solve();

    std::shared_ptr<vecXcd> getQvec() { return qvec; }

    std::shared_ptr<vecXcd> getSols() { return sols; }

    // cmplx getQvec(size_t idx) { return qvec[idx]; }

    void resetSols() { (*sols) = vecXcd::Zero(numSols); }

    void printSols(const std::string&);

private :
    //Solver();
    //Solver(const Solver&) = delete;
    //Solver& operator=(const Solver&) = delete;
    
    std::shared_ptr<Node> root;

    int numSols;
    int maxIter;
    double EPS;

    matXcd Qmat;
    matXcd Hmat;

    vecXcd gvec;
    vecXcd currents;

    double g0;

    std::shared_ptr<vecXcd> qvec;
    std::shared_ptr<vecXcd> sols;
    
};