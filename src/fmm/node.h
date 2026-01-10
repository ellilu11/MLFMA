#pragma once

#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include "../config.h"
#include "../interp.h"
#include "../phys.h"
#include "angles.h"
#include "fmm.h"
#include "tables.h"

class Solver;

using NodeVec = std::vector<std::shared_ptr<FMM::Node>>;
using NodePair = std::pair<std::shared_ptr<FMM::Node>, std::shared_ptr<FMM::Node>>;

class FMM::Node {

public:
    static void initParams(
        const Config&, 
        const std::shared_ptr<Excitation::PlaneWave>&);

    static void linkStates(const std::unique_ptr<Solver>&);

    static void buildTables();

    Node(const SrcVec&, const int, Node* const);

    virtual void initNode() = 0;

    virtual void buildMpoleCoeffs() = 0;

    virtual void buildLocalCoeffs() = 0;

    void printFarSols(const std::string&);

    void printAngles();

    static int getMaxLvl() { return maxLevel; }

    static int getNumNodes() { return numNodes; }

    SrcVec getSrcs() const { return srcs; }
    
    int getBranchIdx() const { return branchIdx; }
    
    double getLeng() const { return nodeLeng; }
    
    int getLevel() const { return level; }

    vec3d getCenter() const { return center; }
    
    Node* getBase() const { return base; }
    
    NodeVec getBranches() const { return branches; }
    
    NodeVec getNbors() const { return nbors; }

    NodeVec getIlist() const { return iList; }

    NodeVec getLeafIlist() const { return leafIlist; }
    
    std::vector<vec2cd> getMpoleCoeffs() const { return coeffs; }
    
    std::vector<vec2cd> getLocalCoeffs() const { return localCoeffs; }

    bool isRoot() const { return base == nullptr; }
    
    template <typename T>
    bool isNodeType() const { return typeid(*this) == typeid(T); }

    bool isSrcless() const { return srcs.empty(); }

    // ========== Test methods ==========
    void testFarfield(double);

    static void printAngularSamples(int);

    static std::shared_ptr<Node> getNode();

    std::vector<vec3cd> getFarSolsFromCoeffs(double);

protected:
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>, const Dir) const;

    void resizeCoeffs();
    
    void buildInteractionList();
    
    void pushSelfToNearNonNbors();

    void buildMpoleToLocalCoeffs();

    void evalLeafIlistSols();
   
    virtual std::shared_ptr<Node> getSelf() = 0;
    
    virtual void buildNeighbors() = 0;

protected:
    static Config config;
    static double wavenum;
    static Angles angles;
    static Tables tables;
    static std::vector<NodePair> nonNearPairs;

    static std::shared_ptr<vecXcd> lvec;
    static std::shared_ptr<vecXcd> rvec;
    static std::shared_ptr<vecXcd> currents;

    inline static int numNodes = 0;
    inline static int maxLevel = 0;

    std::vector<vec2cd> coeffs;
    std::vector<vec2cd> localCoeffs;

    NodeVec branches;
    NodeVec nbors; // list 1
    NodeVec iList; // list 2
    NodeVec leafIlist; // list 4

    SrcVec srcs;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const int level;
    const vec3d center;
};