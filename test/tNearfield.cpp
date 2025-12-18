#include <random>
#include "../src/MLFMA.h"
#include "../src/fileio.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

shared_ptr<Node> Node::getNode() {

    size_t nodeIdx = 0;

    auto node = nodes[nodeIdx];

    // while (node->isNodeType<Stem>() || node->iList.empty() || node->srcs.size() <= 5)
    while (node->isNodeType<Leaf>() || node->isRoot() || node->srcs.size() == 0)
        node = nodes[++nodeIdx];

    cout << "   Selected node " << node->nodeIdx << " at level " << node->level 
         << " of length " << node->nodeLeng
         << " with " << node->iList.size() << " interaction nodes and "
         << node->srcs.size() << " srcs\n\n";

    return node;
}

/*void Node::printAngularSamples(int level) {
    const auto [nth, nph] = getNumAngles(level);

    ofstream thetaFile("out/nf/thetas_nth" + to_string(nth) + ".txt");
    ofstream phiFile("out/nf/phis_nth" + to_string(nth) + ".txt");

    for (int ith = 0; ith < nth; ++ith)
        thetaFile << thetas[level][ith] << '\n';

    for (int iph = 0; iph < nph; ++iph)
        phiFile << phis[level][iph] << '\n';
}*/

void Leaf::testMpoleToLocalInLeaf() {

    auto [nth, nph] = Node::getNumAngles(level);

    // Get sols from local coeffs due to iList (assuming L2L is off)
    ofstream solFile("out/nf/nf_nth"+to_string(nth)+".txt");

    evalFarSols();

    for (const auto& src : srcs)
        solFile << src->getSol() << '\n';

    // Get sols directly from iList
    ofstream solDirFile("out/nf/nf_dir.txt");

    resetSols();

    for (const auto& node : iList)
         evalPairSols(node);

    for (const auto& src : srcs)
        solDirFile << src->getSol() << '\n';

    /* 
    // Print out norms of coeffs
    ofstream coeffFile("out/nf/coeffs.txt");
    for (int idx = 0; idx < nth*nph; ++idx)
        coeffFile << getLocalCoeffs()[idx] << '\n';

    // Print out positions of sources and observers
    ofstream obsFile("out/nf/obss.txt");
    for (const auto& src : srcs)
        obsFile << src->getCenter() << '\n';

    ofstream srcFile("out/nf/srcs.txt");
    for (const auto& src : node->getSrcs())
        srcFile << src->getCenter() << '\n';

    // Print out centers of source and observer box
    std::cout << "   Source box at: " << node->getCenter() << '\n';
    std::cout << "   Target box at: " << center << '\n';
    */
}

void Stem::testShiftedLocalCoeffs() {

    const auto [mth, mph] = getNumAngles(level+1);

    const int branchIdx = 7;

    const auto& inCoeffs = branches[branchIdx]->getMpoleCoeffs();
    const auto& outCoeffs = getShiftedMpoleCoeffs(branchIdx);
    
    ofstream coeffFile("out/coeffs/coeffs.txt");
    ofstream aCoeffFile("out/coeffs/acoeffs.txt");

    for (int idx = 0; idx < mth*mph; ++idx) {
        coeffFile << setprecision(15) << inCoeffs[idx] << '\n';
        aCoeffFile << setprecision(15) << outCoeffs[idx] << '\n';
    }
}

/*void testTransl() {
    // ==================== Test translation ===================== //
    cout << " Testing M2L translations...\n";

    realVec testDists =
        { 2, 2.23607, 2.44949, 2.82843, 3, 3.16228, 3.31662, 3.4641,
          3.60555, 3.74166, 4.12311, 4.24264, 4.3589, 4.69042, 5.19615 };

    const auto& translTable = Node::getTables().transl;

    for (const auto& dist : testDists)
        cout << dist << ' ' << translTable[1].at(dist)[0] << '\n';
}*/

extern auto t = ClockTimes();

int main() {
    // ===================== Read config ==================== //
    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto fmm_start = Clock::now();
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
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

    // ==================== Build tables ===================== //
    cout << " Building angular samples...\n";

    start = Clock::now();

    Node::buildAngularSamples();
    Node::buildTables();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Upward pass ===================== //
    cout << " Computing upward pass...\n";

    start = Clock::now();

    root->buildMpoleCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // L2L test
    auto obsNode = root->getNode();
    auto obsStem = dynamic_pointer_cast<Stem>(obsNode);
    obsStem->testShiftedLocalCoeffs();

    return 0;

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass...\n";
    start = Clock::now();

    root->buildLocalCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Nearfield test =================== //
    // auto obsLevel = obsNode->getLevel();
    // auto [nth, nph] = Node::getNumAngles(obsLevel);
    //ofstream coeffFile("out/nf/lcoeffs_nth"+to_string(nth)+".txt");
    //obsNode->printLocalCoeffs(coeffFile);
    //Node::printAngularSamples(obsLevel);

    // auto obsLeaf = dynamic_pointer_cast<Leaf>(obsNode);
    // obsLeaf->testMpoleToLocalInLeaf();

    return 0;

}