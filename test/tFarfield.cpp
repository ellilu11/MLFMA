#include "../src/MLFMA.h"
#include "../src/interp.h"
#include "../src/leaf.h"
#include "../src/stem.h"

using namespace std;

vec2cd Leaf::getLeafSols(const vec3d R) {
    const int nth = thetas[level].size();
    const int nph = phis[level].size();

    auto rhat = R / R[0];

    vec2cd fld;

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {
        auto weight = thetaWeights[level][ith];

        for (int iph = 0; iph < nph; ++iph) {
            const auto shift = R - center;

            fld += weight * coeffs[idx++] 
                    * Math::expI(wavenum*rhat.dot(shift));
        }
    }

    return fld;
}

void Leaf::testFarFieldFromLeaves(const std::vector<vec3d>& obss) {

    ofstream outFile;
    outFile.open("out/ff.txt");

    for (const auto& leaf : leaves)
        leaf->buildMpoleCoeffs();

    for (const auto& obs : obss) {
        vec2cd fld = vec2cd::Zero();

        for (const auto& leaf : leaves)
            fld = fld + leaf->getLeafSols(obs);

        outFile << fld << '\n';
    }
}

int main() {
    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    auto vertices = importVertices("config/n120/vertices.txt");

    auto tris = importTriangles("config/n120/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    auto srcs = importRWG("config/n120/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    cout << " # Sources:           " << Nsrcs << '\n';
    cout << " Root length:         " << config.rootLeng << '\n';
    cout << " Interpolation order: " << config.order << '\n';
    cout << " Max node RWGs:       " << config.maxNodeSrcs << "\n\n";

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl() << '\n';

    // ==================== Build tables ===================== //
    cout << " Building tables...\n";

    Node::buildAngularSamples();
    Node::buildTables(config);

    // ==================== Test upward pass ===================== //
    cout << " Testing upward pass...\n";

    // Set up observers
    ofstream obsFile;
    obsFile.open("config/obss.txt");

    const int nth = 10;
    const int nph = 20;
    const double R = 50.0 * config.rootLeng;

    std::vector<vec3d> obss;
    for (int ith = 0; ith < nth; ++ith) {
        double th = PI * ith / static_cast<double>(nth);
        for (int iph = 0; iph < nph; ++iph) {
            double ph = 2.0 * PI * iph / static_cast<double>(nph);
            auto obs = vec3d(R, th, ph);
            obss.push_back(obs);
            obsFile << obs << '\n';
        }
    }

    // Print farfield due to leaf coefficients
    Leaf::testFarFieldFromLeaves(obss);

    // root->testFarField(obss);

    // Print direct solution
    ofstream outFileDir;
    outFileDir.open("out/ffDir.txt");

    for (const auto& obs : obss)
        outFileDir << root->getFarSols(obs) << '\n';

    return 0;
}