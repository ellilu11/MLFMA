#pragma once

#include <cassert>
#include <numeric>
#include <queue>
#include "../config.h"
#include "../interp.h"
#include "../phys.h"
#include "angles.h"
#include "fmm.h"
#include "tables.h"

using NodeVec = std::vector<std::shared_ptr<FMM::Node>>;
using NodePair = std::pair<std::shared_ptr<FMM::Node>, std::shared_ptr<FMM::Node>>;

class FMM::Node {
    friend struct Angles;
    friend class Tables;

public:
    static void initStatic(
        const Config&, 
        const std::shared_ptr<Excitation::PlaneWave>&,
        int);

    Node(const SrcVec&, const int, Node* const);

    static void buildTables();

    virtual void buildLists() = 0;

    virtual void resizeCoeffs() = 0;

    virtual void buildMpoleCoeffs() = 0;

    virtual void buildLocalCoeffs() = 0;

    void printFarSols(const std::string&);

    void printAngles();

    static int getMaxLvl() { return maxLevel; }

    static int getNumNodes() { return numNodes; }

    SrcVec getSrcs() const { return srcs; }
    
    int getBranchIdx() const { return branchIdx; }

    Node* getBase() const { return base; }
    
    double getNodeLeng() const { return nodeLeng; }
    
    int getLevel() const { return level; }

    vec3d getCenter() const { return center; }
    
    std::vector<vec2cd> getMpoleCoeffs() const { return coeffs; }
    
    std::vector<vec2cd> getLocalCoeffs() const { return localCoeffs; }

    bool isRoot() const { return base == nullptr; }
    
    bool isSrcless() const { return srcs.empty(); }

    template <typename T>
    bool isNodeType() const { return typeid(*this) == typeid(T); }

    //void testFarfield(double);

    //static std::shared_ptr<Node> getNode();

    //std::vector<vec3cd> getFarSolsFromCoeffs(double);

protected:
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>, const Dir) const;
    
    void buildInteractionList();
    
    void pushSelfToNearNonNbors();

    void buildMpoleToLocalCoeffs();

    void evalLeafIlistSols();
   
    virtual std::shared_ptr<Node> getSelf() = 0;
    
    virtual void buildNeighbors() = 0;

public:
    static std::shared_ptr<vecXcd> lvec;
    static std::shared_ptr<vecXcd> rvec;
    static std::shared_ptr<vecXcd> currents;

protected:
    static Config config;
    static double wavenum;
    static std::vector<Angles> angles;
    static std::vector<Tables> tables;
    static std::vector<NodePair> nonNearPairs;

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