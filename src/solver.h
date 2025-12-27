#pragma once

#include <filesystem>
#include "types.h"
#include "sources/source.h"

class Solver {

public :
    Solver(const SrcVec& srcs);

    void iterateArnoldi(int);

    pair2cd applyGivensRotation(vecXcd&, const vecXcd&, const vecXcd&, int);

    void solve(int numIter);

    cmplx getQvec(size_t idx) { return qvec[idx]; }

    void addToSols(size_t idx, cmplx val) { sols[idx] += val; }

    void resetSols() { sols = vecXcd::Zero(nsols); }

    void printSols(const std::string&);

private :
    //Solver();
    //Solver(const Solver&) = delete;
    //Solver& operator=(const Solver&) = delete;
    
    int nsols;

    matXcd Q;
    matXcd H;

    vecXcd qvec;
    vecXcd sols;
    
    vecXcd gvec;
    vecXcd currents;
};