#include <fstream>
#include <iostream>
#include "../src/MLFMA.h"
#include "../src/interp.h"
#include "../src/leaf.h"
#include "../src/stem.h"

using namespace std;

/* testFarfieldFromLeaves(r)
 * Print total farfield along sampled directions at distance r
 */
void Node::testFarfield(double r) {

    ofstream outFile("out/ff_maxlvl" + to_string(maxLevel) + ".txt");

    outFile << setprecision(15) << scientific;

    const auto [nth, nph] = getNumAngles(level);

    buildMpoleCoeffs();

    auto sols = getFarSolsFromCoeffs(r);

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const vec3d solAbs = sols[idx].cwiseAbs();

            outFile << solAbs << '\n';

            // if (ith == 2) cout << solAbs << '\n';

            idx++;
        }
    }
}

/* testFarfieldFromLeaves(r)
 * Print total farfield along leaf sampled directions at distance r,
 * assuming all non-empty leaves are at same level
 */
void Leaf::testFarfieldFromLeaves(double r) {

    ofstream outFile("out/ff_maxlvl" + to_string(maxLevel) + ".txt");

    outFile << setprecision(15) << scientific;

    const auto [nth, nph] = getNumAngles(maxLevel);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    for (const auto& leaf : leaves) {
        if (leaf->rwgs.empty()) continue;

        // only works if all non-empty leaves are at maxlevel
        assert(leaf->level == maxLevel);

        leaf->buildMpoleCoeffs();

        sols = sols + leaf->getFarSolsFromCoeffs(r);
    }

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const vec3d solAbs = sols[idx].cwiseAbs();

            outFile << solAbs << '\n';

            idx++;
        }
    }
}

/*void Node::testFarfieldDir(double r) {

    ofstream outFile("out/ffDir.txt");

    outFile << setprecision(15) << scientific;

    const auto [nth, nph] = getNumAngles(level);

    auto sols = getFarSols(r);

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const vec3d solAbs = sols[idx].cwiseAbs();

            outFile << solAbs << '\n';

            idx++;
        }
    }
}*/

void Node::printAngularSamples(int level) {
    ofstream thetaFile("out/thetas_lvl" + to_string(level) + ".txt");
    ofstream phiFile("out/phis_lvl" + to_string(level) + ".txt");

    const auto [nth, nph] = getNumAngles(level);

    for (int ith = 0; ith < nth; ++ith)
        thetaFile << thetas[level][ith] << '\n';

    for (int iph = 0; iph < nph; ++iph)
        phiFile << phis[level][iph] << '\n';
}

int main() {
    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    cout << " Constructing RWGs...\n";

    auto vertices = importVertices("config/n120/vertices.txt");

    auto tris = importTriangles("config/n120/faces.txt", vertices, config);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    auto srcs = importRWG("config/n120/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    cout << "   # Sources:       " << Nsrcs << '\n';
    cout << "   Quad precision:  " << static_cast<int>(config.quadPrec) << '\n';
    cout << "   Digit precision: " << config.digits << '\n';
    cout << "   Interp order:    " << config.interpOrder << '\n';
    cout << "   Root length:     " << config.rootLeng << '\n';
    cout << "   Max node RWGs:   " << config.maxNodeSrcs << '\n';
    cout << "   Wave number:     " << Einc->wavenum << '\n';

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

    // ==================== Build tables ===================== //
    cout << "\n Building angular samples...\n";

    Node::buildAngularSamples();
    Node::buildTables();

    // ==================== Test upward pass ===================== //
    const double r = 20.0*config.rootLeng; // pick r >> rootLeng

    // Compute farfield from S2M + M2M
    cout << "\n Computing fields...\n";
    auto start = Clock::now();

    root->testFarfield(r);
    // Leaf::testFarfieldFromLeaves(r);

    auto end = Clock::now();
    Time duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // Print out theta and phi samples at root level
    Node::printAngularSamples(0);

    return 0;
}