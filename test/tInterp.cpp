#include "../src/MLFMA.h"
#include "../src/fileio.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

extern auto t = ClockTimes();

cmplx phiFunc(double ph) {
    return exp(iu*ph);
}

double thetaFunc(double th) {
    return cos(th);
};

cmplx sphFunc(double th, double ph) {
    // return 0.48860251190292 * cos(th);
    // return thetaFunc(th);
    return 0.345494149471336 * sin(th) * exp(iu*ph);
};

void Stem::tInterpPhi(int srcLvl, int tgtLvl) {
    const int order = config.interpOrder;

    // Evaluate function at source nodes
    const auto mph = getNumAngles(srcLvl).second;
    cmplxVec vals;

    for (int iph = 0; iph < mph; ++iph) {
        const double phi = phis[srcLvl][iph];
        vals.push_back(phiFunc(phi));
    }

    // Interpolate over phi
    const auto nph = getNumAngles(tgtLvl).second;
    cmplxVec interpVals(nph, 0.0);
    size_t n = 0;

    for (int jph = 0; jph < nph; ++jph) {
        const auto [interp, nearIdx] = tables.interpPhi[tgtLvl][jph];

        for (int iph = nearIdx+1-order, k = 0; iph <= nearIdx+order; ++iph, ++k) {

            const int iph_wrapped = Math::wrapIdxToRange(iph, mph);

            interpVals[n] += interp[k] * vals[iph_wrapped];
        }

        ++n;
    }

    // Do inner product (weighted)
    const double phiWeight = 2.0*PI / static_cast<double>(nph);
    cmplx intVal = 0.0;
    for (int jph = 0; jph < nph; ++jph) {
        const double phi = phis[tgtLvl][jph];
        intVal += phiWeight * conj(phiFunc(phi)) * interpVals[jph];
    }

    std::cout << "Integrated val using interp: " << std::setprecision(15) << intVal << '\n';
}


void Stem::tAnterpPhi(int srcLvl, int tgtLvl) {
    const int order = config.interpOrder;

    // Evaluate function times weights at source nodes
    const auto mph = getNumAngles(srcLvl).second;
    const auto nph = getNumAngles(tgtLvl).second;
    const int Nph = nph+2*order;
    cmplxVec vals;

    const double phiWeight = 2.0*PI / static_cast<double>(mph);
    for (int iph = 0; iph < mph; ++iph) {
        const double phi = phis[srcLvl][iph];

        vals.push_back(phiWeight*phiFunc(phi));
    }

    // Anterpolate over phi
    cmplxVec extVals(Nph, 0.0);
    for (int iph = 0; iph < mph; ++iph) { // over parent phis anterpolating child phis
        const auto [interp, nearIdx] = tables.interpPhi[srcLvl][iph];

        for (int jph = -order; jph < nph+order; ++jph) { // over child phis to anterpolate

            // shift from jph \in [nearIdx+1-order,nearIdx+order] to k \in [0,2*order-1]
            const int k = jph - (nearIdx+1-order);

            // if iph \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            const int Jph = jph+order;
            extVals[Jph] += interp[k] * vals[iph];
        }
    }

    // Contract extended theta
    cmplxVec anterpVals(nph, 0.0);
    for (int jph = 0; jph < nph; ++jph) {
        const int Jph = jph+order;

        anterpVals[jph] += extVals[Jph];

        if (jph < order || jph >= nph-order) {
            const int Jph_wrapped = Jph + (jph < order ? nph : -nph);
            anterpVals[jph] += extVals[Jph_wrapped];
        }
    }

    // Do inner product
    cmplx intVal = 0.0;
    for (int jph = 0; jph < nph; ++jph) {
        const double phi = phis[tgtLvl][jph];
        intVal += conj(anterpVals[jph]) * phiFunc(phi);
    }

    std::cout << "Integrated val using anterp: " << std::setprecision(15) << intVal << '\n';
}

void Stem::tInterpTheta(int srcLvl, int tgtLvl) {
    const int order = config.interpOrder;

    // Evaluate function at source nodes
    const int mth = getNumAngles(srcLvl).first;
    const int nth = getNumAngles(tgtLvl).first;
    cmplxVec vals;

    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[srcLvl].first[ith];
        vals.push_back(thetaFunc(theta));
    }

    // Interpolate function values to target nodes
    cmplxVec interpVals(nth, 0.0);
 
    for (int jth = 0; jth < nth; ++jth) {

        const auto [interp, nearIdx] = tables.interpTheta[tgtLvl][jth];

        for (int ith = nearIdx+1-order, k = 0; ith <= nearIdx+order; ++ith, ++k) {

            // Flip ith if not in [0, mth-1]
            const int ith_flipped = Math::flipIdxToRange(ith, mth);

            interpVals[jth] += interp[k] * vals[ith_flipped];
        }
    }

    // Do inner product (weighted)
    cmplx intVal = 0.0;
    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[tgtLvl].first[jth];
        const double thetaWeight = thetas[tgtLvl].second[jth];

        intVal += thetaWeight * thetaFunc(theta) * interpVals[jth];
    }

    std::cout << "Integrated val using interp: " << std::setprecision(15) << 2*PI*intVal << '\n';
}

void Stem::tAnterpTheta(int srcLvl, int tgtLvl) {
    const int order = config.interpOrder;

    // Evaluate function times weights at source nodes
    const int mth = getNumAngles(srcLvl).first;
    const int nth = getNumAngles(tgtLvl).first;
    const int Nth = nth+2*order;
    cmplxVec vals;

    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[srcLvl].first[ith];
        const double thetaWeight = thetas[srcLvl].second[ith];
        vals.push_back(thetaWeight*thetaFunc(theta));
    }

    // Anterpolate function values to target nodes on extended grid
    cmplxVec extVals(Nth, 0.0);

    for (int ith = 0; ith < mth; ++ith) { // over parent thetas anterpolating child thetas
        const auto [interp, nearIdx] = tables.interpTheta[srcLvl][ith];

        for (int jth = -order; jth < nth+order; ++jth) { // over child thetas to anterpolate

            // shift from ith \in [t+1-order,t+order] to k \in [0,2*order-1]   
            const int k = jth - (nearIdx+1-order);

            // if ith \notin [t+1-order,t+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            const int Jth = jth+order;
            extVals[Jth] += interp[k] * vals[ith];
        }
    }

    // Contract extended theta
    cmplxVec anterpVals(nth, 0.0);
    for (int jth = 0; jth < nth; ++jth) {
        const int Jth = jth+order;

        anterpVals[jth] += extVals[Jth];

        // Handle nodes near poles
        if (jth < order || jth >= nth-order) {
            const int jth_flipped = (jth < order ? -jth-1 : 2*nth-jth-1);

            anterpVals[jth] += extVals[jth_flipped+order];

            // std::cout << std::setprecision(9) << extVals[jth_flipped+order] << '\n';
        }
        //
    }

    // Do inner product
    cmplx intVal = 0.0;
    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[tgtLvl].first[jth];
        intVal += anterpVals[jth] * thetaFunc(theta);
        // std::cout << std::setprecision(6) << anterpVals[jth] << ' ';
    }
    // std::cout << '\n';

    std::cout << "Integrated val using anterp: " << std::setprecision(15) << 2.0*PI*intVal << '\n';
}

void Stem::tInterp(int srcLvl, int tgtLvl) {
    const int order = config.interpOrder;

    // Evaluate function at source nodes
    const auto [mth, mph] = getNumAngles(srcLvl);
    const auto [nth, nph] = getNumAngles(tgtLvl);
    cmplxVec vals;

    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[srcLvl].first[ith];

        for (int iph = 0; iph < mph; ++iph) {
            const double phi = phis[srcLvl][iph];

            vals.push_back(sphFunc(theta, phi));
        }
    }

    // Interpolate function values to target nodes
    cmplxVec interpVals(nth*nph, 0.0);
    addInterpCoeffs(vals, interpVals, srcLvl, tgtLvl);

    // Do inner product (weighted)
    const double phiWeight = 2.0*PI / static_cast<double>(nph);
    cmplx intVal = 0.0;
    size_t l = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[tgtLvl].first[jth];
        const double thetaWeight = thetas[tgtLvl].second[jth];

        for (int jph = 0; jph < nph; ++jph) {
            const double phi = phis[tgtLvl][jph];

            // intVal += thetaWeight * phiWeight * conj(sphFunc(theta, phi)) * interpVals[l];

            const cmplx prod = thetaWeight * phiWeight * conj(sphFunc(theta, phi)) * interpVals[l];
            intVal += prod;

            // std::cout << std::setprecision(15)
            //    << jth << ' ' << jph << ' ' << theta << ' ' << phi << ' ' << prod << '\n';

            // std::cout << std::setprecision(15) 
            //    << sphFunc(theta, phi) << ' ' << interpVals[l] << ' ' << prod << '\n';

            ++l;
        }
    }
    assert(l == interpVals.size());

    std::cout << "Integrated val using interp: " << std::setprecision(15) << intVal << '\n';
}

void Stem::tAnterp(int srcLvl, int tgtLvl) {
    const int order = config.interpOrder;

    // Evaluate function times weights at source nodes
    const auto [mth, mph] = getNumAngles(srcLvl);
    const auto [nth, nph] = getNumAngles(tgtLvl);
    cmplxVec vals;

    const double phiWeight = 2.0*PI / static_cast<double>(mph);
    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[srcLvl].first[ith];
        const double thetaWeight = thetas[srcLvl].second[ith];

        for (int iph = 0; iph < mph; ++iph) {
            const double phi = phis[srcLvl][iph];

            vals.push_back(phiWeight*thetaWeight*sphFunc(theta, phi));
        }
    }

    // Anterpolate function values to target nodes on extended grid
    cmplxVec anterpVals(nth*nph, 0.0);
    addAnterpCoeffs(vals, anterpVals, srcLvl, tgtLvl);

    // Do inner product
    cmplx intVal = 0.0;
    size_t l = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[tgtLvl].first[jth];

        for (int jph = 0; jph < nph; ++jph) {
            const double phi = phis[tgtLvl][jph];

            const cmplx prod = conj(anterpVals[l]) * sphFunc(theta, phi);
            intVal += prod;

            //std::cout << std::setprecision(15)
            //    << sphFunc(theta, phi) << ' ' << anterpVals[l] << ' ' << prod << '\n';

            // std::cout << std::setprecision(15) << jth << ' ' << jph << ' ' << theta << ' ' << phi << ' ' << prod << '\n';

            ++l;
        }
    }

    assert(l == anterpVals.size());

    std::cout << "Integrated val using anterp: " << std::setprecision(15) << intVal << '\n';
}

int main() {
    // ===================== Read config ==================== //
    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    Node::initParams(config, Einc);

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto fmm_start = Clock::now();
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

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
    root->initNode();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // Do test
    const int lvl = 2;
    // Stem::tInterpPhi(lvl+1,lvl);
    // Stem::tAnterpPhi(lvl,lvl+1);
    // Stem::tInterpTheta(lvl+1,lvl);
    // Stem::tAnterpTheta(lvl,lvl+1);
    Stem::tInterp(lvl+1,lvl);
    Stem::tAnterp(lvl,lvl+1);

    return 0;
}

/*void Stem::testInvInterp(int level) {

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

                const bool nearPole = ith != ith_flipped;

                int iph_shifted = iph;

                if (nearPole)
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

}*/