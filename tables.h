#pragma once

#include "math.h"

struct Tables {
    Tables() = default;
    Tables(const int order, const Precision prec) {
        buildYlmTables(order);
        buildQuadTables(prec);
        buildExpTables(order);
    }

    void buildYlmTables(const int);
    void buildDirTables(const Precision);

    // Ylm tables
    std::vector<realVec> coeffYlm_;
    std::vector<realVec> fallingFact_;
    std::vector<realVec> legendreSum_;
    std::vector<realVec> fracCoeffYlm_;
    std::vector<realVec> A_;
    std::vector<realVec> Aexp_;

    // quad tables
    std::vector<pair2d> quadCoeffs_;
    std::vector<int> quadLengs_;
};

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
}

void Tables::buildDirTables(const Precision prec) {

}
