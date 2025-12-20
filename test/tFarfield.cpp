#include <fstream>
#include <iostream>
#include "../src/MLFMA.h"
#include "../src/fileio.h"
#include "../src/interp.h"
#include "../src/leaf.h"
#include "../src/stem.h"

using namespace std;

/* getFarSolsFromCoeffs(r)
 * Return sols at all sampled directions at distance r
 * from mpole coefficients of this node
 */
std::vector<vec3cd> Node::getFarSolsFromCoeffs(double r) {
    assert(!srcs.empty());

    const cmplx kernel = C * wavenum * exp(iu*wavenum*r) / r;

    const auto [nth, nph] = getNumAngles(level);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto kvec = tables.khat[level][idx] * wavenum;

            sols[idx] =
                kernel * exp(-iu*kvec.dot(center))
                // * tables.matFromSph[level][idx]
                * coeffs[idx];

            idx++;
        }
    }

    return sols;
}

/* testFarfield(r)
 * Print total farfield along sampled directions at distance r
 */
void Node::testFarfield(double r) {

    const auto [nthRoot, nphRoot] = getNumAngles(0);

    ofstream outFile("out/ff/ff_maxlvl" + to_string(maxLevel) + "_nth" + to_string(nthRoot) + ".txt");

    outFile << setprecision(15) << scientific;

    const auto [nth, nph] = getNumAngles(level);

    buildMpoleCoeffs();

    auto sols = getFarSolsFromCoeffs(r);

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const vec3d solAbs = sols[idx++].cwiseAbs();

            outFile << solAbs << '\n';

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

/* testFarfieldFromLeaves(r)
    * Print total farfield along leaf sampled directions at distance r,
    * assuming all non-empty leaves are at same level
    */
    /*void TestLeaf::testFarfieldFromLeaves(double r) {

            ofstream outFile("out/ff/ff_maxlvl" + to_string(maxLevel) + ".txt");

            outFile << setprecision(15) << scientific;

            const auto [nth, nph] = getNumAngles(maxLevel);

            std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

            for (const auto& leaf : leaves) {
                if (leaf->srcs.empty()) continue;

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
        }*/


void Node::printAngularSamples(int level) {

    const auto [nthRoot, nphRoot] = getNumAngles(0);

    ofstream thetaFile("out/ff/thetas_nth" + to_string(nthRoot) + ".txt");
    ofstream phiFile("out/ff/phis_nph" + to_string(nthRoot) + ".txt");

    const auto [nth, nph] = getNumAngles(level);

    for (int ith = 0; ith < nth; ++ith)
        thetaFile << thetas[level][ith] << '\n';

    for (int iph = 0; iph < nph; ++iph)
        phiFile << phis[level][iph] << '\n';
}

extern auto t = ClockTimes();

int main() {
    // ===================== Read config ==================== //
    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);

    auto nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
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
    const double r = 10.0*config.rootLeng; // pick r >> rootLeng

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