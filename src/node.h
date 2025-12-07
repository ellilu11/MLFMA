#pragma once

#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include "clock.h"
#include "config.h"
#include "interp.h"
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
    // static const int getExpansionOrder() { return order; }

    // static void setExpansionOrder(const int p) { order = p; }
    
    static int getPrecision() { return prec; }

    static int getMaxLvl() { return maxLevel; }
    
    static int getNumNodes() { return numNodes; }

    static pair2i getNumAngles(const int level) { 
        return std::make_pair(thetas[level].size(), phis[level].size());
    }

    static void setNodeParams(const Config&, const std::shared_ptr<Src>&);

    static void buildAngularSamples();

    static void buildTables(const Config&);

public:
    RWGVec getRWGs() const { return rwgs; }
    
    const int getBranchIdx() const { return branchIdx; }
    
    const double getLeng() const { return nodeLeng; }
    
    const vec3d getCenter() const { return center; }
    
    Node* getBase() const { return base; }
    
    const NodeVec getBranches() const { return branches; }
    
    NodeVec getNbors() const { return nbors; }

    NodeVec getIlist() const { return iList; }

    NodeVec getLeafIlist() const { return leafIlist; }
    
    std::vector<vec2cd> getMpoleCoeffs() const { return coeffs; }
    
    std::vector<vec2cd> getLocalCoeffs() const { return localCoeffs; }

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

    void buildMpoleToLocalCoeffs();
   
    std::vector<vec2cd> getShiftedLocalCoeffs(const int) const;

    void evalLeafIlistSols();

    void evalPairSols(const std::shared_ptr<Node>&);

    void evalSelfSols();

    vec2cd getFarSols(const vec3d);
   
    virtual std::shared_ptr<Node> getSelf() = 0;
    
    virtual void buildNeighbors() = 0;

    virtual void buildLists() = 0;
    
    virtual void buildMpoleCoeffs() = 0;
    
    virtual void buildLocalCoeffs() = 0;
    
    virtual void printNode(std::ofstream&) = 0;

    // ========== Test methods ==========
    // void testFarField(const int, const int);

protected:
    static int order;
    static int prec;
    static int maxNodeSrcs;
    static int maxLevel;
    static double rootLeng;
    static double wavenum;
    static std::vector<realVec> phis;
    static std::vector<realVec> thetas;
    static std::vector<realVec> thetaWeights;
    static std::vector<int> Ls;
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
    NodeVec iList; // list 2
    NodeVec leafIlist; // list 4

    std::vector<vec2cd> coeffs;
    std::vector<vec2cd> localCoeffs;

    int label;
};