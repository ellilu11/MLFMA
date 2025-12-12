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

// TODO: move into config
constexpr double c0 = 299792458.0; 
constexpr double mu0 = 1.256637E-6; 
constexpr double q = 3.5; // TODO: pick an optimal value

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

    static void setNodeParams(const Config&, const std::shared_ptr<Src>&);

    static void buildAngularSamples();

    static void buildTables() { tables = Tables(0); }

public:
    RWGVec getRWGs() const { return rwgs; }
    
    int getBranchIdx() const { return branchIdx; }
    
    double getLeng() const { return nodeLeng; }
    
    vec3d getCenter() const { return center; }
    
    Node* getBase() const { return base; }
    
    NodeVec getBranches() const { return branches; }
    
    NodeVec getNbors() const { return nbors; }

    NodeVec getIlist() const { return iList; }

    NodeVec getLeafIlist() const { return leafIlist; }
    
    std::vector<vec3cd> getMpoleCoeffs() const { return coeffs; }

    std::pair<vec3cd, vec3cd> getPolarCoeffs() const { return polarCoeffs; }
    
    std::vector<vec3cd> getLocalCoeffs() const { return localCoeffs; }

    bool isRoot() const { return base == nullptr; }
    
    template <typename T>
    bool isNodeType() const { return typeid(*this) == typeid(T); }

    bool isSrcless() const { return rwgs.empty(); }

    // void resetSols() { for (const auto& p : rwg) p->resetSol(); }

public:
    Node(const RWGVec&, const int, Node* const);
    
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>&, const Dir) const;
    
    void buildInteractionList();
    
    void pushSelfToNearNonNbors();

    void buildMpoleToLocalCoeffs();

    void evalLeafIlistSols();

    void evalPairSols(const std::shared_ptr<Node>&);

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
    // std::vector<vec3d> getObssAtAngularSamples(double);

    void testFarfield(double);

    void testFarfieldDir(double);

    static void printAngularSamples(int);

    static Tables getTables(){ return tables; }

protected:
    static Config config;
    static double wavenum;
    static int maxLevel;
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

    // TODO: template the vec type
    std::vector<vec3cd> coeffs; 
    std::pair<vec3cd, vec3cd> polarCoeffs;
    std::vector<vec3cd> localCoeffs;

    int label;
};