#include "farfield.h"

FMM::Farfield::Farfield(const std::shared_ptr<FMM::Node>& root) {
    if (root->isLeaf()) return;

    buildLevels();
    buildGlRadPats();
    resizeCoeffs(root);
}

/* buildLevels()
 * Build angular samples and interpolation/translation tables for each level
 */
void FMM::Farfield::buildLevels() {
    std::cout << " Building FMM operators...        ";
    auto start = Clock::now();

    Level::buildDists();

    levels.reserve(maxLevel+1);
    for (int lvl = 0; lvl <= maxLevel; ++lvl)
        levels.emplace_back(lvl);
    
    Level::clearDists();

    for (int lvl = 0; lvl < maxLevel; ++lvl)
        levels[lvl].buildInterpTables(levels[lvl+1]);

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
}

/* buildGlRadPats()
 * Build plane wave expansions due to sources in leaf nodes
 */
void FMM::Farfield::buildGlRadPats() {
    std::cout << " Building plane wave expansions...";
    auto start = Clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (int iLeaf = 0; iLeaf < glLeaves.size(); ++iLeaf)
        buildRadPats(glLeaves[iLeaf]);

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
}

/* resizeCoeffs(node)
 * Resize mpole and local coeffs of node and its branches to match 
 * number of directions at each level
 */
void FMM::Farfield::resizeCoeffs(const std::shared_ptr<FMM::Node>& node) {
    size_t nDir = levels[node->lvl].getNumDirs();
    node->coeffs.resize(nDir);
    node->localCoeffs.resize(nDir);

    for (const auto& branch : node->branches)
        resizeCoeffs(branch);
}

/* buildRadPats()
 * Build radiation patterns due to sources in leaf node
 */
void FMM::Farfield::buildRadPats(const std::shared_ptr<FMM::Node>& node) {
    if (node->isSrcless()) return;
    const Level& level = levels[node->lvl];
    size_t nDir = level.getNumDirs();

    size_t nSrc = node->srcs.size(), iSrc = 0;
    node->radPats.resize(nSrc);
    if (config.ie != IE::EFIE) node->recPatsH.resize(nSrc); 

    for (const auto& src : node->srcs) {
        Coeffs radPat(nDir);
        Coeffs recPatH(nDir); // TODO: Don't resize to nDir if EFIE

        for (int iDir = 0; iDir < nDir; ++iDir) {
            vec3d kvec = level.khat[iDir] * config.k;
            const mat23d& toThPh = level.toThPh[iDir];

            auto [rad, radNormal] = src->getRadsAlongDir(node->center, kvec);
            radPat.setCoeffAlongDir(toThPh * rad, iDir);

            if (config.ie == IE::EFIE) continue;

            recPatH.setCoeffAlongDir(
                toThPh * iu * (kvec.cross(radNormal).conjugate()),
                iDir);
        }

        node->radPats[iSrc] = std::move(radPat);
        if (config.ie != IE::EFIE)
            node->recPatsH[iSrc] = std::move(recPatH);

        ++iSrc;
    }
}

/* buildMpoleCoeffs()
 * (S2M) Build multipole coefficients from sources in this leaf
 */
void FMM::Farfield::buildMpoleCoeffs(const std::shared_ptr<FMM::Node>& node) {
    Coeffs& coeffs = node->coeffs;

    // Zero out coeffs before accumulation
    coeffs.fillZero();
    if (node->isSrcless() || node->isRoot()) return;

    size_t nDir = levels[node->lvl].getNumDirs();
    Eigen::Map<arrXcd> coeffsTheta(coeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> coeffsPhi(coeffs.phi.data(), nDir);

    // #pragma omp parallel for reduction(+:coeffsTheta,coeffsPhi)
    for (int iSrc = 0; iSrc < node->srcs.size(); ++iSrc) {
        const auto& src = node->srcs[iSrc];
        Coeffs& radPat = node->radPats[iSrc];

        Eigen::Map<arrXcd> radPatTheta(radPat.theta.data(), nDir);
        Eigen::Map<arrXcd> radPatPhi(radPat.phi.data(), nDir);

        coeffsTheta += Solver::lvec[src->getIdx()] * radPatTheta;
        coeffsPhi += Solver::lvec[src->getIdx()] * radPatPhi;
    }
}

/* buildMpoleCoeffs(node, isStem = true)
 * (M2M) Build mpole coeffs by merging branch mpole coeffs
 */
void FMM::Farfield::buildMpoleCoeffs(const std::shared_ptr<FMM::Node>& node, bool isStem) {
    if (node->isLeaf()) return;

    int order = config.interpOrder;
    int lvl = node->lvl;
    size_t mDir = levels[lvl+1].getNumDirs();

    // Zero out coeffs before merging
    node->coeffs.fillZero();

    for (const auto& branch : node->branches) {
        if (branch->isSrcless()) continue;

        // Shift branch coeffs to center of this node
        vec3d dX = node->center - branch->center;

        Coeffs shiftedCoeffs(mDir);
        for (int iDir = 0; iDir < mDir; ++iDir) {
            const vec3d& kvec = levels[lvl+1].khat[iDir] * config.k;
            cmplx shift = exp(iu*kvec.dot(dX));

            shiftedCoeffs.theta[iDir] = shift * branch->coeffs.theta[iDir];
            shiftedCoeffs.phi[iDir] = shift * branch->coeffs.phi[iDir];
        }

        // Interpolate shifted coeffs to this node's angular grid
        addInterpCoeffs(shiftedCoeffs, node->coeffs, lvl+1, lvl);
    }
}

/* translateCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center,
 * then apply integration weights for anterpolation
 */
void FMM::Farfield::translateCoeffs(const std::shared_ptr<FMM::Node>& node) {
    if (node->iList.empty()) return;

    Coeffs& localCoeffs = node->localCoeffs;
    localCoeffs.fillZero();
    size_t nDir = localCoeffs.size();

    // Translate mpole coeffs into local coeffs
    const VecHashMap<arrXcd>& transl = levels[node->lvl].transl;
    Eigen::Map<arrXcd> localTheta(localCoeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> localPhi(localCoeffs.phi.data(), nDir);

    for (const auto& srcNode : node->iList) {
        assert(!srcNode->isSrcless());

        vec3d dX = node->center - srcNode->center;
        arrXcd transl_dX = transl.at(dX/node->nodeLeng);

        Eigen::Map<arrXcd> mpoleTheta(srcNode->coeffs.theta.data(), nDir);
        Eigen::Map<arrXcd> mpolePhi(srcNode->coeffs.phi.data(), nDir);

        localTheta += transl_dX * mpoleTheta;
        localPhi += transl_dX * mpolePhi;
    }

    // Apply integration weights
    // TODO: Absorb weights into translation table
    const Level& level = levels[node->lvl];
    auto [nth, nph] = level.getNumAngles();
    double phiWeight = 2.0*PI / static_cast<double>(nph);

    size_t iDir = 0;
    for (int ith = 0; ith < nth; ++ith) {
        double thetaWeight = level.weights[ith];

        for (int iph = 0; iph < nph; ++iph) {
            localCoeffs.theta[iDir] *= thetaWeight * phiWeight;
            localCoeffs.phi[iDir] *= thetaWeight * phiWeight;
            ++iDir;
        }
    }
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
FMM::Coeffs FMM::Farfield::getShiftedLocalCoeffs(
    FMM::Node* node, int branchIdx) const 
{
    int lvl = node->lvl;
    size_t mDir = levels[lvl].getNumDirs();
    size_t nDir = levels[lvl+1].getNumDirs();

    Coeffs outCoeffs(nDir);
    if (node->iList.empty()) return outCoeffs;

    // Shift local coeffs to center of branch
    vec3d dX = node->branches[branchIdx]->center - node->center;

    Coeffs shiftedCoeffs(mDir);
    for (int iDir = 0; iDir < mDir; ++iDir) {
        vec3d kvec = levels[lvl].khat[iDir] * config.k;
        cmplx shift = exp(iu*kvec.dot(dX));

        shiftedCoeffs.theta[iDir] = shift * node->localCoeffs.theta[iDir];
        shiftedCoeffs.phi[iDir] = shift * node->localCoeffs.phi[iDir];
    }

    // Anterpolate shifted coeffs to branch's angular grid
    addAnterpCoeffs(shiftedCoeffs, outCoeffs, lvl, lvl+1);

    return outCoeffs;
}

/* buildLocalCoeffs()
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void FMM::Farfield::buildLocalCoeffs(const std::shared_ptr<FMM::Node>& node) {
    if (node->isSrcless()) return;

    if (!node->isRoot()) {
        // Add L2L to M2L
        if (!node->base->isRoot()) // TODO: Avoid adding Coeffs directly
            node->localCoeffs += 
                getShiftedLocalCoeffs(node->base, node->branchIdx);
    }
}

/* evalFarSols()
 * (L2T) Evaluate sols from local expansion due to far nodes
 */
void FMM::Farfield::evalFarSols(const std::shared_ptr<FMM::Node>& node) {
    int lvl = node->lvl;
    if (node->isSrcless() || lvl <= 1) return;

    size_t nDir = levels[lvl].getNumDirs();
    Eigen::Map<arrXcd> localTheta(node->localCoeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> localPhi(node->localCoeffs.phi.data(), nDir);

    auto getIntRad = [&](Coeffs& rec) {
        Eigen::Map<arrXcd> recTheta(rec.theta.data(), nDir);
        Eigen::Map<arrXcd> recPhi(rec.phi.data(), nDir);
        return (recTheta.conjugate() * localTheta +
            recPhi.conjugate() * localPhi).sum();
    };

    size_t iObs = 0;
    for (const auto& obs : node->srcs) {
        switch (config.ie) {
            case IE::EFIE: {
                cmplx intRadE = getIntRad(node->radPats[iObs]);
                Solver::rvec[obs->getIdx()] += 4.0 * PI * config.C_efie * intRadE;
                break;
            }
            case IE::MFIE: {
                cmplx intRadH = getIntRad(node->recPatsH[iObs]);
                Solver::rvec[obs->getIdx()] += 4.0 * PI * config.C_mfie * intRadH;
                break;
            }
            case IE::CFIE: {
                cmplx intRadE = getIntRad(node->radPats[iObs]);
                cmplx intRadH = getIntRad(node->recPatsH[iObs]);
                Solver::rvec[obs->getIdx()] += 
                    4.0 * PI * (config.C_efie * intRadE + config.C_mfie * intRadH);
                break;
            }
        }
        ++iObs;
    }
}

/* evaluateSols()
 * Do FMM main loop to evaluate farfield contributions to sols, and add to rvec
 */
void FMM::Farfield::evaluateSols() {
    // S2M
    auto start = Clock::now();
    #pragma omp parallel for schedule(dynamic)
    for (int iLeaf = 0; iLeaf < glLeaves.size(); ++iLeaf)
        buildMpoleCoeffs(glLeaves[iLeaf]);
    t.S2M += Clock::now() - start;

    // M2M
    start = Clock::now();
    for (int lvl = maxLevel-1; lvl >= 0; --lvl) {
        auto& nodes = glNodesByLvl[lvl];

        #pragma omp parallel for if (lvl >= config.pivotLvl)
        for (int iNode = 0; iNode < nodes.size(); ++iNode)
            buildMpoleCoeffs(nodes[iNode], true);
    }
    t.M2M += Clock::now() - start;

    // M2L
    start = Clock::now(); 
    #pragma omp parallel for schedule(dynamic, 1)
    for (int iNode = 0; iNode < glNodes.size(); ++iNode)
        translateCoeffs(glNodes[iNode]);
    t.M2L += Clock::now() - start;

    // L2L
    start = Clock::now();
    for (int lvl = 1; lvl <= maxLevel; ++lvl) {
        auto& nodes = glNodesByLvl[lvl];

        #pragma omp parallel for if (lvl >= config.pivotLvl)
        for (int iNode = 0; iNode < nodes.size(); ++iNode)
            buildLocalCoeffs(nodes[iNode]);
    }
    t.L2L += Clock::now() - start;

    // L2T
    start = Clock::now();
    #pragma omp parallel for schedule(dynamic)
    for (int iLeaf = 0; iLeaf < glLeaves.size(); ++iLeaf)
        evalFarSols(glLeaves[iLeaf]);
    t.L2T += Clock::now() - start;
}