#pragma once

#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include "clock.h"
#include "config.h"
#include "interp.h"
#include "tables.h"
#include "sources/rwg.h"

extern ClockTimes t;

constexpr int DIM = 3;
constexpr int numDir = 26;

// TODO: move into config
constexpr double c0 = 299792458.0; 
constexpr double mu0 = 1.256637E-6;
constexpr double Q = 3.5; // TODO: pick an optimal value

const cmplx C = -iu * c0 * mu0 / (4.0 * PI); // TODO: Find a better place to put this

enum class Dir {
    W, E, S, N, D, U,
    SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
    DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
};

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {
    friend struct Tables;

public:
    static int getMaxLvl() { return maxLevel; }
    
    static int getNumNodes() { return numNodes; }

    static pair2i getNumAngles(const int level) { 
        return std::make_pair(thetas[level].size(), phis[level].size());
    }

    static void setNodeParams(const Config&, const std::shared_ptr<PlaneWave>&);

    static void buildAngularSamples();

    static void buildTables() { tables = Tables(config); }

public:
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
    
    std::vector<vec3cd> getMpoleCoeffs() const { return coeffs; }

    // std::vector<vec3cd> getRadPats() const { return radPats; } //

    std::pair<vec3cd, vec3cd> getPolarCoeffs() const { return polarCoeffs; }
    
    std::vector<vec3cd> getLocalCoeffs() const { return localCoeffs; }

    bool isRoot() const { return base == nullptr; }
    
    template <typename T>
    bool isNodeType() const { return typeid(*this) == typeid(T); }

    bool isSrcless() const { return srcs.empty(); }

    void resetSols() { for (const auto& src : srcs) src->resetSol(); }

public:
    Node(const SrcVec&, const int, Node* const);
    
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>, const Dir) const;
    
    void buildInteractionList();
    
    void pushSelfToNearNonNbors();

    void buildMpoleToLocalCoeffs();

    void evalLeafIlistSols();

    void evalPairSols(const std::shared_ptr<Node>);

    void evalSelfSols();

    std::vector<vec3cd> getFarSolsFromCoeffs(double);

    std::vector<vec3cd> getFarSols(double);
   
    virtual std::shared_ptr<Node> getSelf() = 0;
    
    virtual void buildNeighbors() = 0;

    virtual void buildLists() = 0;
    
    virtual void buildMpoleCoeffs() = 0;
    
    virtual void buildLocalCoeffs() = 0;
    
    virtual void printNode(std::ofstream&) = 0;

    // ========== Test methods ==========
    static Tables getTables() { return tables; }

    void testFarfield(double);

    // void testFarfieldDir(double);

    static void printAngularSamples(int);

    static std::shared_ptr<Node> getNode();

    virtual void printLocalCoeffs(std::ofstream& f) = 0;

protected:
    static Config config;
    static double wavenum;
    inline static int numNodes = 0;
    inline static int maxLevel = 0;

    static std::vector<realVec> thetas;
    static std::vector<realVec> thetaWeights;
    static std::vector<realVec> phis;
    static std::vector<int> Ls;
    
    static Tables tables;

    // TODO: template the vec type
    std::vector<vec3cd> coeffs;
    std::pair<vec3cd, vec3cd> polarCoeffs;
    std::vector<vec3cd> localCoeffs;

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

    // Test members
    static NodeVec nodes;
    int nodeIdx;
    int label;
};