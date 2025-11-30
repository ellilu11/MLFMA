#pragma once

#include "math.h"

struct Tables {
    Tables() = default;
    Tables(const int maxLevel, const double k, std::vector<realVec>& thetas) {
        buildDirTables(maxLevel, k, thetas);
        // buildYlmTables(order);
    }


    void buildDirTables(const int, const double, std::vector<realVec>&);
    // void buildYlmTables(const int);

    // dir tables
    std::vector<std::vector<mat3d>> ImRR;
    std::vector<std::vector<vec3d>> kvec;

    // Ylm tables
    /*std::vector<realVec> coeffYlm_;
    std::vector<realVec> fallingFact_;
    std::vector<realVec> legendreSum_;
    std::vector<realVec> fracCoeffYlm_;
    std::vector<realVec> A_;
    std::vector<realVec> Aexp_;*/
};

void Tables::buildDirTables(
    const int maxLevel, const double k, std::vector<realVec>& thetas) {
    for (int level = 0; level <= maxLevel; ++level) {

        const int nth = thetas[level].size();
        const int nph = 2*nth;

        std::vector<mat3d> ImRR_lvl;
        std::vector<vec3d> kvec_lvl;

        for (int ith = 0; ith < nth; ++ith) {
            const double th = thetas[level][ith];
            for (int iph = 0; iph < nph; ++iph) {
                const double ph = 2.0*PI*iph/static_cast<double>(nph);

                ImRR_lvl.push_back(Math::IminusRR(th, ph));
                kvec_lvl.push_back(Math::vecSph(k, th, ph));
            }
        }

        ImRR.push_back(ImRR_lvl);
        kvec.push_back(kvec_lvl);
    }
}

/*
void Tables::buildYlmTables(const int order) {
    auto binom = [](double x, int k) {
        return Math::fallingFactorial(x, k) / Math::factorial(k);
        };

    for (int l = 0; l <= 2*order; ++l) {
        realVec coeffYlm_l, fallingFact_l, legendreSum_l, fracCoeffYlm_l, A_l, Aexp_l;

        for (int m = 0; m <= l; ++m) {
            coeffYlm_l.push_back(Math::coeffYlm(l, m));
            fallingFact_l.push_back(Math::fallingFactorial(l, m));
            legendreSum_l.push_back(binom(l, m) * binom((l+m-1)/2.0, l));
            fracCoeffYlm_l.push_back(sqrt((l-m)/static_cast<double>(l+m)));
        }

        auto pm_l = Math::pm(l);
        for (int m = -l; m <= l; ++m) {
            A_l.push_back(pm_l /
                std::sqrt(static_cast<double>(Math::factorial(l-m)*Math::factorial(l+m))));
            Aexp_l.push_back(1.0 /
                std::sqrt(static_cast<double>(Math::factorial(l-m)*Math::factorial(l+m))));
        }

        coeffYlm_.push_back(coeffYlm_l);
        fallingFact_.push_back(fallingFact_l);
        legendreSum_.push_back(legendreSum_l);
        fracCoeffYlm_.push_back(fracCoeffYlm_l);
        A_.push_back(A_l);
        Aexp_.push_back(Aexp_l);
    }
}*/
