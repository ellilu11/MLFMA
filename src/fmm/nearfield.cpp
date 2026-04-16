#include "nearfield.h"
#include "../mesh/tripairs.h"

FMM::Nearfield::Nearfield(size_t nsrcs)
    : nearMat(nsrcs, nsrcs)
{
    findNodePairs();
    buildTriPairs();

    buildNearMatrix();
    Mesh::glPairsToIdx = {};
    Mesh::glTriPairs.clear();
}

/* findNodePairs()
 * From list of leaves, find all near neighbor leaf pairs
 */
void FMM::Nearfield::findNodePairs() {
    for (const auto& leaf : glLeaves) {
        selfPairs.emplace_back(leaf, leaf);

        for (const auto& nbor : leaf->nearNbors) {
            assert(nbor->isLeaf());
            if (leaf < nbor) nearPairs.emplace_back(leaf, nbor);
        }
    }

    for (const auto& pair : glNonNearPairs) {
        const auto& [obsNode, srcNode] = pair;
        nearPairs.emplace_back(obsNode, srcNode);
    }
}

/* findTriPairs()
 * From list of node pairs, find all near neighbor triangle pairs
 * and populate Mesh::glTriPairs
 */
void FMM::Nearfield::buildTriPairs() {
    std::cout << " Building triangle moments...     ";
    auto start = Clock::now();

    size_t iPair = 0;

    // Assign indices to each self pair
    for (const auto& [leaf, srcLeaf] : selfPairs) {
        assert(leaf == srcLeaf);
        const auto& iTris = leaf->iTris;

        for (auto iTri0 : iTris)
            for (auto iTri1 : iTris) {
                if (iTri0 > iTri1) continue;

                pair2i pair(iTri0, iTri1);
                Mesh::glPairsToIdx.emplace(pair, iPair++);
            }
    }

    // Assign indices to each near pair
    for (const auto& [obsLeaf, srcNode] : nearPairs) {
        const auto& iTris0 = obsLeaf->iTris, & iTris1 = srcNode->iTris;

        for (auto iTri0 : iTris0)
            for (auto iTri1 : iTris1) {
                pair2i pair = std::minmax(iTri0, iTri1);
                Mesh::glPairsToIdx.emplace(pair, iPair++);
            }
    }

    // Construct triangle moments of all found pairs
    Mesh::glTriPairs = Mesh::TriPairs(iPair);

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
}

/* getNearCapacity()
 * Get number of nonzero entries in near matrix to reserve space for triplets
 */
size_t FMM::Nearfield::getNearCapacity() {
    return
        std::accumulate(nearPairs.begin(), nearPairs.end(), 0,
            [](size_t sum, const auto& nearPair) {
                const auto& [obsLeaf, srcNode] = nearPair;
                size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();
                return sum + 2*nObss*nSrcs;
            }
        ) +
        std::accumulate(selfPairs.begin(), selfPairs.end(), 0,
            [](size_t sum, const auto& selfPair) {
                const auto& [leaf, srcLeaf] = selfPair;
                size_t nSrcs = leaf->srcs.size();
                return sum + nSrcs*nSrcs;
            }
        );
}

/* buildNearMatrix()
 * From list of node pairs, build near matrix by computing pairwise contributions
 * between sources in each node pair and adding to nearMat
 */
//
void FMM::Nearfield::buildNearMatrix() {
    std::cout << " Building nearfield matrix...     ";
    auto start = Clock::now();

    std::vector<Eigen::Triplet<cmplx>> trips;
    trips.reserve(getNearCapacity());

    // Build pair-node contributions to near matrix
    std::vector<std::vector<Eigen::Triplet<cmplx>>> tripsPair(config.numThreads);
    #pragma omp parallel num_threads(config.numThreads) 
    {
        int tid = omp_get_thread_num();
        auto& local = tripsPair[tid];

        #pragma omp for
        for (int iPair = 0; iPair < nearPairs.size(); ++iPair) {
            const auto& [obsLeaf, srcNode] = nearPairs[iPair];
            size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();

            for (int iObs = 0; iObs < nObss; ++iObs) {
                auto obs = obsLeaf->srcs[iObs];
                size_t obsIdx = obs->getIdx(); // global index of obs

                for (int iSrc = 0; iSrc < nSrcs; ++iSrc) {
                    auto src = srcNode->srcs[iSrc];
                    size_t srcIdx = src->getIdx(); // global index of src

                    double mass = obs->getIntegratedMass(src);
                    cmplx efie = config.C_efie * obs->getIntegratedEFIE(src),
                        mfieObs = config.C_mfie * (obs->getIntegratedMFIE(src) + mass),
                        mfieSrc = config.C_mfie * (src->getIntegratedMFIE(obs) + mass);

                    local.emplace_back(obsIdx, srcIdx, efie+mfieObs);
                    local.emplace_back(srcIdx, obsIdx, efie+mfieSrc);
                }
            }
        }
    }
    for (const auto& local : tripsPair)
        trips.insert(trips.end(), local.begin(), local.end());

    // Build self-node contributions to near matrix
    std::vector<std::vector<Eigen::Triplet<cmplx>>> tripsSelf(config.numThreads);
    #pragma omp parallel num_threads(config.numThreads)
    {
        int tid = omp_get_thread_num();
        auto& local = tripsSelf[tid];

        #pragma omp for
        for (int iPair = 0; iPair < selfPairs.size(); ++iPair) {
            const auto& [leaf, srcLeaf] = selfPairs[iPair];
            assert(leaf == srcLeaf);
            size_t nSrcs = leaf->srcs.size();

            for (int iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
                auto obs = leaf->srcs[iObs];
                size_t obsIdx = obs->getIdx(); // global index of obs

                for (int iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                    auto src = leaf->srcs[iSrc];
                    size_t srcIdx = src->getIdx(); // global index of src

                    double mass = obs->getIntegratedMass(src);
                    cmplx efie = config.C_efie * obs->getIntegratedEFIE(src);
                    cmplx mfieObs = config.C_mfie * (obs->getIntegratedMFIE(src) + mass);

                    local.emplace_back(obsIdx, srcIdx, efie+mfieObs);

                    if (iSrc != iObs) { // Only add self-term contribution once!
                        cmplx mfieSrc = config.C_mfie * (src->getIntegratedMFIE(obs) + mass);
                        local.emplace_back(srcIdx, obsIdx, efie+mfieSrc);
                    }
                }
            }
        }
    }
    for (const auto& local : tripsSelf)
        trips.insert(trips.end(), local.begin(), local.end());

    nearMat.setFromTriplets(trips.begin(), trips.end());
    nearMat.makeCompressed();

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
}

/* evaluateSols()
 * (S2T) Multiply near matrix by lvec and add to rvec to get nearfield contribution to rvec
 */
void FMM::Nearfield::evaluateSols() {
    auto start = Clock::now();

    assert(nearMat.cols() == Solver::lvec.rows());
    assert(nearMat.rows() == Solver::rvec.rows());

    Solver::rvec += nearMat * Solver::lvec;

    t.S2T += Clock::now() - start;
}

/* printNearMatrix()
 * Print near matrix to file for debugging
 */
void FMM::Nearfield::printNearMatrix(const std::string& fname) const {
    std::ofstream of(fname);
    matXcd denseMat(nearMat);
    for (int i = 0; i < denseMat.rows(); ++i) {
        for (int j = 0; j < denseMat.cols(); ++j)
            of << denseMat(i, j).real() << ' ';
        of << '\n';
    }
}