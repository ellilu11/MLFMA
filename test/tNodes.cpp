#include <fstream>
#include <iostream>
#include "../src/MLFMA.h"
#include "../src/interp.h"
#include "../src/leaf.h"
#include "../src/stem.h"



int main() {
    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    cout << " Constructing RWGs...\n";

    auto vertices = importVertices("config/nlarge/vertices.txt");

    auto tris = importTriangles("config/nlarge/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    auto srcs = importRWG("config/nlarge/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    cout << "   # Sources:     " << Nsrcs << '\n';
    cout << "   Root length:   " << config.rootLeng << '\n';
    cout << "   Interp order:  " << config.order << '\n';
    cout << "   Max node RWGs: " << config.maxNodeSrcs << "\n";
    cout << "   Wave number:   " << Einc->wavenum << "\n\n";

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    const int maxLevel = Node::getMaxLvl();

    root->buildLists();

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << maxLevel << '\n';

    // ==================== Print nodes ===================== //

    ofstream file("out/nodes.txt");

    Leaf::printLeaves(file);

    return 0;
}