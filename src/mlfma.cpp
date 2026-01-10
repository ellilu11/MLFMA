#include <fstream>
#include "MLFMA.h"
#include "clock.h"
#include "config.h"
#include "fileio.h"
#include "fmm/fmm.h"

using namespace std;

extern auto t = ClockTimes();

int main() {
    // ===================== Read config ==================== //
    cout << " Importing sources...\n";

    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    FMM::Node::initParams(config, Einc);

    // ==================== Set up nodes ==================== //
    cout << " Setting up nodes...\n";

    shared_ptr<FMM::Node> root;
    if (nsrcs > config.maxNodeSrcs)
        root = make_shared<FMM::Stem>(srcs, 0, nullptr);
    else
        root = make_shared<FMM::Leaf>(srcs, 0, nullptr);

    cout << "   # Nodes: " << FMM::Node::getNumNodes() << '\n';
    // cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << FMM::Node::getMaxLvl() << "\n\n";

    // ==================== Build tables ===================== //
    cout << " Building tables...\n";

    auto start = Clock::now();

    FMM::Node::buildTables();
    root->initNode();

    auto end = Clock::now();
    Time duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build nearfield ===================== //
    cout << " Building nearfield interactions...\n";

    start = Clock::now();

    FMM::Leaf::buildNearRads();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build radpats ===================== //
    cout << " Building radiation patterns...\n";

    start = Clock::now();

    FMM::Leaf::buildRadPats();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Solve iterative FMM ================ //
    cout << " Solving w/ FMM...\n";

    constexpr int MAX_ITER = 1000;
    constexpr double EPS = 1.0E-6;

    auto solver = make_unique<Solver>(srcs, root, MAX_ITER, EPS);
    FMM::Node::linkStates(solver);

    start = Clock::now();

    solver->evalRvec(0);

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Total elapsed time: " << duration_ms.count() << " ms\n\n";

    solver->printSols("rvec.txt");
    //root->printFarSols("ff_nq7.txt");
    //root->printAngles();

    if (!config.evalDirect) return 0;

    // ================== Solve iterative direct ================ //
    FMM::Leaf::resetLeaves();
    root = make_shared<FMM::Leaf>(srcs, 0, nullptr);
    root->initNode();

    cout << " Building nearfield interactions...\n";

    start = Clock::now();
    FMM::Leaf::buildNearRads();
    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    cout << " Solving w/ direct...\n";
    solver = make_unique<Solver>(srcs, root, MAX_ITER, EPS);
    FMM::Node::linkStates(solver);

    solver->evalRvec(0);

    solver->printSols("rvecDir.txt");

    return 0;
}