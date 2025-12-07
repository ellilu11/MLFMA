#include <fstream>
#include <iostream>
#include "clock.h"
#include "config.h"
#include "MLFMA.h"

using namespace std;

extern auto t = ClockTimes();

int main() {
    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    auto vertices = importVertices("config/n120/vertices.txt");

    auto tris = importTriangles("config/n120/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    auto srcs = importRWG("config/n120/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config,Einc);

    cout << " # Sources:           " << Nsrcs << '\n';
    cout << " Root length:         " << config.rootLeng << '\n';
    cout << " Interpolation order: " << config.order << '\n';
    // cout << " Exponential order:   " << Node::getExponentialOrder() << '\n';
    cout << " Max node RWGs:       " << config.maxNodeSrcs << "\n\n";

    // return 0;

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto fmm_start = Clock::now();
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    auto end = Clock::now();
    Time duration_ms = end - start;

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl() << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";
    
    // return 0;

    // ==================== Build tables ===================== //
    cout << " Building tables...\n";

    start = Clock::now();

    Node::buildAngularSamples();
    Node::buildTables(config);

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";
    
    // return 0;

    // ==================== Upward pass ===================== //
    cout << " Computing upward pass...\n";

    start = Clock::now();

    root->buildMpoleCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    return 0;

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass...\n";
    start = Clock::now();

    root->buildLocalCoeffs();

    end = Clock::now();
    duration_ms = end - start;

    // ================== Evaluate solutions ================= //
    cout << " Evaluating solutions...\n";
    start = Clock::now();

    Leaf::evaluateSols();

    end = Clock::now();
    duration_ms = end - start;
    Time fmm_duration_ms = end - fmm_start;

    cout << " FMM total elapsed time: " << fmm_duration_ms.count() << " ms\n";

}