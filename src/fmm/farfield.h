#pragma once

#include "node.h"

class FMM::Farfield {

public:
    Farfield(const std::shared_ptr<Node>&);

    void evaluateSols();

private:
    void buildLevels();

    void buildGlRadPats();

    void buildRadPats(const std::shared_ptr<Node>&);

    void resizeCoeffs(const std::shared_ptr<Node>&);

    void buildMpoleCoeffs(const std::shared_ptr<Node>&);

    void buildMpoleCoeffs(const std::shared_ptr<Node>&, bool);

    void translateCoeffs(const std::shared_ptr<Node>&);

    Coeffs getShiftedLocalCoeffs(Node*, int) const;

    void buildLocalCoeffs(const std::shared_ptr<Node>&);

    void evalFarSols(const std::shared_ptr<Node>&);

    void addInterpCoeffs(const Coeffs&, Coeffs&, int, int) const;

    void addAnterpCoeffs(const Coeffs&, Coeffs&, int, int) const;

    void addInterpCoeffAlongTh(const NodeVec&, int, int, int) const;

private:
    std::vector<Level> levels;
};