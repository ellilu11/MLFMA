#include <fstream>
#include <iostream>
#include "clock.h"
#include "config.h"
#include "mlfma.h"

using namespace std;

int main() {

    // ===== Tests ===== //
    //for (int l = 1; l <= 5; ++l) {
    //    auto [nodes, weights] = Math::gaussLegendre(l, 1.0E-9, 0.0, PI);
    //    cout << "l = " << l << " ";
    //    for (int k = 0; k < l; ++k)
    //        cout << '(' << nodes[k] << ',' << weights[k] << ") ";
    //    cout << '\n';
    //}
    //return 0;

    // ==================== Import geometry ==================== //
    Config config("config/config.txt");

    vector<vec3d> vertices = importVertices("config/vertices.txt");

    // cout << vertices.size() << '\n';

    TriVec tris = importTriangles("config/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    RWGVec srcs = importRWG("config/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config,Einc);

    // cout << " Source file:       " << fpath.generic_string() << '\n';
    cout << " # Sources:         " << Nsrcs << '\n';
    cout << " Root length:       " << config.rootLeng << '\n';
    cout << " Expansion order:   " << config.order << '\n';
    cout << " Exponential order: " << Node::getExponentialOrder() << '\n';
    cout << " Max node RWGs:    "  << config.maxNodeSrcs << '\n' << '\n';

    return 0;

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl << '\n' << '\n';

    // ============ Construct directional quantities ========= //

    Node::buildThetaSamples();
    Node::buildTables(config);

    return 0;

    // ==================== Upward pass ===================== //
    cout << " Computing upward pass...\n";

    root->buildMpoleCoeffs();

    return 0;
}