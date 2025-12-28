#pragma once

#include <filesystem>
#include "node.h"
#include "types.h"
#include "sources/source.h"

//namespace Solver {
//    class Solver;
//    struct StateVecs;
//};

class Solver {

public :
    Solver(const SrcVec& srcs);

    void iterateArnoldi(int);

    pair2cd applyGivensRotation(vecXcd&, const vecXcd&, const vecXcd&, int);

    void updateCurrent(int, int);

    void solve(int);

    std::shared_ptr<vecXcd> getQvec() { return qvec; }

    std::shared_ptr<vecXcd> getSols() { return sols; }

    // cmplx getQvec(size_t idx) { return qvec[idx]; }

    // void addToSols(size_t idx, cmplx val) { sols[idx] += val; }

    void resetSols() { (*sols) = vecXcd::Zero(nsols); }

    void printSols(const std::string&);

private :
    //Solver();
    //Solver(const Solver&) = delete;
    //Solver& operator=(const Solver&) = delete;
    
    int nsols;

    matXcd Q;
    matXcd H;

    std::shared_ptr<vecXcd> qvec;
    std::shared_ptr<vecXcd> sols;
    
    vecXcd gvec;
    vecXcd currents;
};