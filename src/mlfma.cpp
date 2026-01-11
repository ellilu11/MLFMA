#include <fstream>
#include "MLFMA.h"
#include "clock.h"
#include "config.h"
#include "fileio.h"
#include "fmm/fmm.h"

using namespace FMM;

extern auto t = ClockTimes();

int main() {
    // ===================== Read config ==================== //
    std::cout << " Importing sources...\n";

    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    Node::initStatic(config, Einc, nsrcs);

    // ==================== Set up nodes ==================== //
    std::cout << " Setting up nodes...\n";

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    std::cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    std::cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    std::cout << "   Max node level: " << Node::getMaxLvl() << "\n\n";

    // ==================== Build nearfield ===================== //
    std::cout << " Building nearfield interactions...\n";

    auto start = Clock::now();
    Leaf::buildNearRads();
    auto end = Clock::now();
    Time duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build tables ======================== //
    std::cout << " Building translation operators...\n";

    start = Clock::now();
    Node::buildTables();
    root->resizeCoeffs(); // TODO: Hide this call
    end = Clock::now();
    duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build radpats ===================== //
    std::cout << " Building radiation patterns...\n";

    start = Clock::now();
    Leaf::buildRadPats();
    end = Clock::now();
    duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Solve iterative FMM ================ //
    std::cout << " Solving w/ FMM...\n";

    constexpr int MAX_ITER = 500;
    constexpr double EPS = 1.0E-6;

    auto solver = make_unique<Solver>(srcs, root, MAX_ITER, EPS,
        Node::lvec, Node::rvec, Node::currents);

    start = Clock::now();
    solver->solve();
    end = Clock::now();
    duration_ms = end - start;
    std::cout << "   Total elapsed time: " << duration_ms.count() << " ms\n\n";

    solver->printSols("curr_nq7.txt");
    //root->printFarSols("ff_nq7.txt");

    if (!config.evalDirect) return 0;

    // ================== Solve iterative direct ================ //
    Node::initStatic(config, Einc, nsrcs);

    root = make_shared<Leaf>(srcs, 0, nullptr);
    root->buildLists();

    std::cout << " Building nearfield interactions...\n";

    start = Clock::now();
    Leaf::buildNearRads();
    end = Clock::now();
    
    duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    std::cout << " Solving w/ direct...\n";
    solver = make_unique<Solver>(srcs, root, MAX_ITER, EPS,
        Node::lvec, Node::rvec, Node::currents);

    solver->solve();

    solver->printSols("currDir_nq7.txt");

    return 0;
}