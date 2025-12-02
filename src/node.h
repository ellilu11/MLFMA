#pragma once

#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include "clock.h"
#include "config.h"
#include "rwg.h"
#include "tables.h"

constexpr int DIM = 3;
constexpr int numDir = 26;

enum class Dir {
    W, E, S, N, D, U,
    SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
    DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
};

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {

public:
    static const int getExpansionOrder() { return order; }

    static void setExpansionOrder(const int p) { order = p; }
    
    static const int getExponentialOrder() { return orderExp; }

    static const int getMaxLvl() { return maxLevel; }
    
    static const int getNumNodes() { return numNodes; }

    static void setNodeParams(const Config&, const std::shared_ptr<Src>&);

    static void buildThetaSamples();

    static void buildTables(const Config&);

public:
    RWGVec getRWGs() const { return rwgs; }
    
    const int getBranchIdx() const { return branchIdx; }
    
    const double getLeng() const { return nodeLeng; }
    
    const vec3d getCenter() const { return center; }
    
    Node* getBase() const { return base; }
    
    const NodeVec getBranches() const { return branches; }
    
    NodeVec getNbors() const { return nbors; }

    NodeVec getDirList(const int dir) const { return dirList[dir]; }

    NodeVec getLeafIlist() const { return leafIlist; }
    
    std::vector<vec3cd> getMpoleCoeffs() const { return coeffs; }
    
    std::vector<vec3cd> getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }
    
    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    // void resetSols() { for (const auto& p : rwg) p->resetSol(); }

public:
    Node(const RWGVec&, const int, Node* const);
    
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>&, const Dir) const;
    
    void buildInteractionList();
    
    void pushSelfToNearNonNbors();
   
    std::vector<vec2cd> getShiftedLocalCoeffs(const int) const;

    void evalLeafIlistSols();

    void evalPairSols(const std::shared_ptr<Node>&);

    void evalSelfSols();
   
    virtual std::shared_ptr<Node> getSelf() = 0;
    
    virtual void buildNeighbors() = 0;

    virtual void buildLists() = 0;
    
    virtual void buildMpoleCoeffs() = 0;
    
    virtual void propagateExpCoeffs() = 0;
    
    virtual void buildLocalCoeffs() = 0;
    
    virtual void printNode(std::ofstream&) = 0;

protected:
    static int order;
    static int orderExp;
    static int maxNodeSrcs;
    static int maxLevel;
    static double rootLeng;
    static double wavenum;
    static std::vector<realVec> phis;
    static std::vector<realVec> thetas;
    static std::vector<realVec> thetaWeights;
    static Tables tables;
    inline static int numNodes = 0;

    RWGVec rwgs;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const int level;
    const vec3d center;

    NodeVec branches;
    NodeVec nbors; // list 1
    std::array<NodeVec,6> dirList; // list 2, indexed by direction
    NodeVec leafIlist; // list 4

    std::vector<vec3cd> coeffs;
    std::vector<vec3cd> localCoeffs;

    int label;
};