#pragma once

#include <filesystem>
#include "types.h"
#include "sources/source.h"

class Node;

class Solver {

public :
    Solver(SrcVec& srcs, std::shared_ptr<Node>, int, double);

    void evalRvec(int);

    void iterateArnoldi(int);

    void applyGivensRotation(vecXcd&, vecXcd&, vecXcd&, int);

    void updateGvec(vecXcd&, vecXcd&, int);

    void solve();

    std::shared_ptr<vecXcd> getLvec() { return lvec; }

    std::shared_ptr<vecXcd> getRvec() { return rvec; }

    std::shared_ptr<vecXcd> getCurrents() { return currents; }

    void resetRvec() { (*rvec) = vecXcd::Zero(numSrcs); }

    void printSols(const std::string&);

private :
    //Solver();
    //Solver(const Solver&) = delete;
    //Solver& operator=(const Solver&) = delete;
    
    std::shared_ptr<Node> root;

    int numSrcs;
    int maxIter;
    double EPS;

    matXcd Qmat;
    matXcd Hmat;
    vecXcd gvec;
    double g0;

    std::shared_ptr<vecXcd> lvec;
    std::shared_ptr<vecXcd> rvec;
    std::shared_ptr<vecXcd> currents;
  
};