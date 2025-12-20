#include "../src/MLFMA.h"
#include "../src/fileio.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

extern auto t = ClockTimes();

void Stem::testInvInterp(int level) {

    auto sphFunc = [](double th, double ph) {
        return sin(th) * exp(iu*ph);
    };

    // Evaluate function at source nodes
    const auto [mth, mph] = getNumAngles(level);
    cmplxVec vals;

    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[level][ith];

        for (int iph = 0; iph < mph; ++iph) {
            const double phi = phis[level][iph];

            vals.push_back(sphFunc(theta, phi));
        }
    }

    // Interpolate function values to target nodes
    const auto [nth, nph] = getNumAngles(level+1);
    const int order = config.interpOrder;
    cmplxVec innerVals(nth*mph, 0.0);

    size_t m = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const int t = tables.invIdxTheta[level][jth];

        for (int iph = 0; iph < mph; ++iph) {

            for (int ith = t+1-order, k = 0; ith <= t+order; ++ith, ++k) {

                const int ith_flipped = Math::flipIdxToRange(ith, mth);

                const bool outOfRange = ith != ith_flipped;

                int iph_shifted = iph;

                if (outOfRange)
                    iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                const int m_shifted = ith_flipped*mph + iph_shifted;

                innerVals[m] +=
                    tables.invInterpTheta[level][jth][k] * vals[m_shifted];
            }

            m++;
        }
    }

    cmplxVec interpedVals(nth*nph, 0.0);
    size_t n = 0;
    for (int jth = 0; jth < nth; ++jth) {

        for (int jph = 0; jph < nph; ++jph) {
            const int s = tables.invIdxPhi[level][jph]; 

            for (int iph = s+1-order, k = 0; iph <= s+order; ++iph, ++k) {

                const int iph_wrapped = Math::wrapIdxToRange(iph, mph);

                interpedVals[n] +=
                    tables.invInterpPhi[level][jph][k]
                    * innerVals[jth*mph + iph_wrapped];
            }

            n++;
        }
    }

    ofstream outFile0("out/valsInterped.txt");
    for (int idx = 0; idx < nth*nph; ++idx)
        outFile0 << interpedVals[idx].real() << ' ' << interpedVals[idx].imag() << '\n';

    // Evaluate function at target nodes directly
    ofstream outFile1("out/valsDirect.txt");

    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[level+1][jth];

        for (int jph = 0; jph < nph; ++jph) {
            const double phi = phis[level+1][jph];

            const cmplx val = sphFunc(theta, phi);

            outFile1 << val.real() << ' ' << val.imag() << '\n';
        }
    }

}

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

    // Do test
    Stem::testInvInterp(0);

    return 0;
}