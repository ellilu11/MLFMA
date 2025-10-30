#include <fstream>
#include <iostream>
#include "clock.h"
#include "config.h"
#include "mlfma.h"

using namespace std;

int main() {
    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    vector<vec3d> vertices = importVertices("config/vertices.txt");

    cout << vertices.size() << '\n';

    TriVec tris = importTriangles("config/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    RWGVec srcs = importRWG("config/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    return 0;

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

}