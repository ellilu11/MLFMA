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
    auto start = Clock::now();

    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    auto solver = make_shared<Solver>(srcs);

    Node::initNodes(config, Einc, solver);

    auto end = Clock::now();
    Time duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Set up domain ==================== //
    cout << " Setting up nodes...\n";
    auto fmm_start = Clock::now();
    start = Clock::now();

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    end = Clock::now();
    duration_ms = end - start;

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl() << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build tables ===================== //
    cout << " Building tables...\n";

    start = Clock::now();

    Node::buildAngularSamples();
    Node::buildTables();

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

    // ==================== Upward pass ===================== //
    cout << " Computing upward pass...\n";

    start = Clock::now();

    root->buildMpoleCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (S2M): " << t.S2M.count() << " ms\n";
    cout << "   Elapsed time (M2M): " << t.M2M.count() << " ms\n\n";

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass...\n";
    start = Clock::now();

    root->buildLocalCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (M2L): " << t.M2L.count() << " ms\n";
    cout << "   Elapsed time (L2L): " << t.L2L.count() << " ms\n\n";

    // ================== Evaluate solutions ================= //
    cout << " Evaluating solutions...\n";
    start = Clock::now();

    Leaf::evaluateSols();

    end = Clock::now();
    duration_ms = end - start;
    Time fmm_duration_ms = end - fmm_start;

    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (L2T): " << t.L2T.count() << " ms\n";
    cout << "   Elapsed time (S2T): " << t.S2T.count() << " ms\n\n";
    cout << " FMM total elapsed time: " << fmm_duration_ms.count() << " ms\n";

    // solver->printSols("sol_d" + to_string(config.digits) + ".txt");
    solver->printSols("sol.txt");

    if (!config.evalDirect) return 0;

    // ================== Compute direct ===================== //
    solver->resetSols();

    cout << "\n Computing direct...\n";
    start = Clock::now();

    root->evalSelfSols();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    solver->printSols("solDir.txt");

    return 0;
}