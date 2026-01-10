#pragma once

#include "fmm.h"

struct FMM::Angles {
    Angles() = default;

    Angles(double, double, int, int);

    pair2i getNumAngles(const int level) const {
        return std::make_pair(thetas[level].size(), phis[level].size());
    }

    size_t getNumAllAngles(const int level) const {
        return thetas[level].size() * phis[level].size();
    }

    std::vector<realVec> thetas;
    std::vector<realVec> thetaWeights;
    std::vector<realVec> phis;
    std::vector<int> Ls;
};

FMM::Angles::Angles(double wavenum, double rootLeng, int digits, int maxLevel) {

    thetas.resize(maxLevel+1);
    thetaWeights.resize(maxLevel+1);
    phis.resize(maxLevel+1);
    Ls.resize(maxLevel+1);

    std::cout << "   (Lvl,Nth,Nph) =\n";

    for (int level = 0; level <= maxLevel; ++level) {
        const double nodeLeng = rootLeng / pow(2.0, level);

        // Use excess bandwidth formula
        const int tau = ceil((1.73*wavenum*nodeLeng +
            2.16*pow(digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));

        Ls[level] = floor(0.50*tau); // TODO: Find optimal formula

        // Construct thetas
        const int nth = tau+1;
        auto [nodes, weights] = Interp::gaussLegendre(nth, EPS_NR, 0.0, PI);

        // Absorb sin(theta) into weights
        std::transform(weights.begin(), weights.end(), nodes.begin(), weights.begin(),
            [](double weight, double theta) { return weight * sin(theta); }
        );

        thetas[level] = nodes;
        thetaWeights[level] = weights;

        // Construct phis
        const int nph = 2*nth;
        realVec phis_lvl(nph);

        for (int iph = 0; iph < nph; ++iph)
            phis_lvl[iph] = 2.0*PI*iph/static_cast<double>(nph);

        phis[level] = phis_lvl;

        std::cout << "   (" << level << "," << nth << "," << nph << ")\n";
    }
}