#include <fstream>
#include <iostream>
#include "clock.h"
#include "config.h"
#include "mlfma.h"

using namespace std;

int main() {

    // ===== Tests ===== //
    for (int l = 1; l <= 5; ++l) {
        auto [nodes, weights] = Math::gaussLegendreTheta(l, 1.0E-9);
        cout << "l = " << l << " ";
        for (int k = 0; k < l; ++k)
            cout << '(' << nodes[k] << ',' << weights[k] << ") ";
        cout << '\n';
    }

    return 0;

    // ==================== Import geometry ==================== //
    Config config("config/config.txt");

    vector<vec3d> vertices = importVertices("config/vertices.txt");

    cout << vertices.size() << '\n';

    TriVec tris = importTriangles("config/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    RWGVec srcs = importRWG("config/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config);

    return 0;

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    Node::setThetaSamples(config);

}