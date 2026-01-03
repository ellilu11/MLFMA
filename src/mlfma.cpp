#include <fstream>
#include "MLFMA.h"
#include "clock.h"
#include "config.h"
#include "fileio.h"

using namespace std;

extern auto t = ClockTimes();

int main() {
    // ===================== Read config ==================== //
    cout << " Importing sources...\n";

    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    Node::initParams(config, Einc);

    // ==================== Set up nodes ==================== //
    cout << " Setting up nodes...\n";

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    // cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl() << "\n\n";

    // ==================== Build tables ===================== //
    cout << " Building tables...\n";

    auto start = Clock::now();

    Node::buildAngularSamples();
    Node::buildTables();
    root->initNode();

    auto end = Clock::now();
    Time duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build nearfield ===================== //
    cout << " Building nearfield interactions...\n";

    start = Clock::now();

    Leaf::buildNearRads();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build radpats ===================== //
    cout << " Building radiation patterns...\n";

    start = Clock::now();

    Leaf::buildRadPats();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Solve iterative FMM ================ //
    cout << " Solving w/ FMM...\n";

    constexpr int MAX_ITER = 1000;
    constexpr double EPS = 1.0E-6;

    auto solver = make_unique<Solver>(srcs, root, MAX_ITER, EPS);
    Node::linkStates(solver);

    start = Clock::now();

    solver->solve();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Total elapsed time: " << duration_ms.count() << " ms\n\n";

    solver->printSols("curr_nq7.txt");
    root->printFarSols("ff_nq7.txt");
    root->printAngles();

    if (!config.evalDirect) return 0;

    // ================== Solve iterative direct ================ //
    Leaf::resetLeaves();
    root = make_shared<Leaf>(srcs, 0, nullptr);
    root->initNode();

    cout << " Building nearfield interactions...\n";

    start = Clock::now();
    Leaf::buildNearRads();
    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    cout << " Solving w/ direct...\n";
    auto solverDir = make_unique<Solver>(srcs, root, MAX_ITER, EPS);
    Node::linkStates(solverDir);

    start = Clock::now();
    solverDir->solve();
    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    solverDir->printSols("currDir_nq7.txt");

    return 0;
}