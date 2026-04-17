#include "tripairs.h"

Mesh::TriPairs::TriPairs(size_t nPair) :
    pairsToIdx(glPairsToIdx.begin(), glPairsToIdx.end()),
    nPair(nPair)
{
    assert(nPair > 0);
    buildNumCommon();

    if (config.ie == IE::EFIE || config.ie == IE::CFIE)
        buildMomentsEFIE();

    if (config.ie == IE::MFIE || config.ie == IE::CFIE)
        buildMomentsMFIE();

    buildIntegratedInvR();

    if (config.ie == IE::MFIE || config.ie == IE::CFIE)
        buildIntegratedInvRcubed();
}

/* buildNumCommon()
 * Compute number of common vertices between all triangle pairs
 */
void Mesh::TriPairs::buildNumCommon() {
    nCommons.resize(nPair);

    for (const auto& [iTris, iPair] : glPairsToIdx) {
        if (iTris.first == iTris.second) {
            nCommons[iPair] = 3;
            continue;
        }

        const auto [tri0, tri1] = getTriPair(iTris);
        int nCommon = 0;

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (tri0.iVerts[i] == tri1.iVerts[j])
                    ++nCommon;
        assert(nCommon < 3);
        nCommons[iPair] = nCommon;
    }
}

/* buildMomentsEFIE()
 * Build EFIE moments for all triangle pairs
 */
void Mesh::TriPairs::buildMomentsEFIE() {
    momentsEFIE.resize(nPair);
    double k = config.k;

    #pragma omp parallel for num_threads(config.numThreads)
    for (int i = 0; i < pairsToIdx.size(); ++i) {
        const auto& [iTris, iPair] = pairsToIdx[i];

        auto& [m00, m10, m01, m11] = momentsEFIE[iPair];
        int nCommon = nCommons[iPair];

        const auto& [obsTri, srcTri] = getTriPair(iTris);

        for (const auto& [obs, obsWeight] : obsTri.quads) {
            for (const auto& [src, srcWeight] : srcTri.quads) {
                double r = (obs-src).norm();

                cmplx G = obsWeight*srcWeight / (4.0*PI); // apply 1/(4pi) factor
                if (nCommon >= nCommonThres)
                    G *= (lzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r);
                else {
                    assert(!lzero(r));
                    G *= exp(iu*k*r) / r;
                }

                m00 += G;
                m10 += obs * G;
                m01 += src * G;
                m11 += obs.dot(src) * G;
            }
        }
    }
}

/* buildMomentsMFIE()
 * Build N-MFIE moments for all triangle pairs
 */
void Mesh::TriPairs::buildMomentsMFIE() {
    momentsMFIE.resize(nPair);
    momentsMFIE2.resize(nPair);
    double k = config.k, k2 = k*k;

    #pragma omp parallel for num_threads(config.numThreads)
    for (int i = 0; i < pairsToIdx.size(); ++i) {
        const auto& [iTris, iPair] = pairsToIdx[i];

        auto& [m000, m001, m10, m01, m11] = momentsMFIE[iPair];
        auto& [n000, n001, n10, n01, n11] = momentsMFIE2[iPair];
        int nCommon = nCommons[iPair];

        const auto& [obsTri, srcTri] = getTriPair(iTris);

        vec3d obsNhat = obsTri.nhat, srcNhat = srcTri.nhat;

        for (const auto& [obs, obsWeight] : obsTri.quads) {
            // For common triangles, use analytic integration of -1/2 J term
            if (nCommon == 3) continue;

            for (const auto& [src, srcWeight] : srcTri.quads) {
                const vec3d& R = obs-src;
                double r = R.norm(), r2 = r*r, r3 = r*r2;
                assert(!lzero(r));

                vec3cd gradG = obsWeight*srcWeight * R / (4.0*PI*r3); // apply 1/(4pi) factor

                if (nCommon >= nCommonThres)
                    gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r) + 0.5*k2*r2 + 1.0); // double check signs
                else
                    gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);

                // Overall minus sign from flipping J x gradG to gradG x J
                m000 -= -obsNhat.dot(gradG);
                m001 -= gradG;
                m10 -= (obs.dot(gradG) * obsNhat - obsNhat.dot(gradG) * obs);
                m01 -= obsNhat.cross(gradG.cross(src));
                m11 -= obs.dot(obsNhat.cross(gradG.cross(src)));

                n000 += -srcNhat.dot(gradG);
                n001 += gradG;
                n10 += (src.dot(gradG) * srcNhat - srcNhat.dot(gradG) * src);
                n01 += srcNhat.cross(gradG.cross(obs));
                n11 += src.dot(srcNhat.cross(gradG.cross(obs)));
            }
        }
    }
}

/* buildIntegratedInvR()
 * Build integrated 1/R and its symmetric case for triangle pair
 */
void Mesh::TriPairs::buildIntegratedInvR() {
    intsInvR.resize(nPair);
    intsInvR2.resize(nPair);

    #pragma omp parallel for num_threads(config.numThreads)
    for (int i = 0; i < pairsToIdx.size(); ++i) {
        const auto& [iTris, iPair] = pairsToIdx[i];
        if (nCommons[iPair] < nCommonThres) continue;

        const auto [obsTri, srcTri] = getTriPair(iTris);

        intRads& intInvR = intsInvR[iPair];
        intInvR.reserve(Triangle::numQuads);
        for (const auto& [obs, weight] : obsTri.quads)
            intInvR.emplace_back(srcTri.getIntegratedInvR(obs));

        intRads& intInvR2 = intsInvR2[iPair];
        intInvR2.reserve(Triangle::numQuads);
        for (const auto& [src, weight] : srcTri.quads)
            intInvR2.emplace_back(obsTri.getIntegratedInvR(src));
    }
}

/* buildIntegratedInvRcubed()
 * Build integrated 1/R^3 and its symmetric case for triangle pair
 */
void Mesh::TriPairs::buildIntegratedInvRcubed() {
    intsInvRcubed.resize(nPair);
    intsInvRcubed2.resize(nPair);

    #pragma omp parallel for num_threads(config.numThreads)
    for (int i = 0; i < pairsToIdx.size(); ++i) {
        const auto& [iTris, iPair] = pairsToIdx[i];
        if (nCommons[iPair] < nCommonThres) continue;

        const auto [obsTri, srcTri] = getTriPair(iTris);

        intRads& intInvR = intsInvRcubed[iPair];
        intInvR.reserve(Triangle::numQuads);
        for (const auto& [obs, weight] : obsTri.quads)
            intInvR.emplace_back(srcTri.getIntegratedInvRcubed(obs));

        intRads& intInvR2 = intsInvRcubed2[iPair];
        intInvR2.reserve(Triangle::numQuads);
        for (const auto& [src, weight] : srcTri.quads)
            intInvR2.emplace_back(obsTri.getIntegratedInvRcubed(src));
    }
}

/*
void Mesh::TriPairs::buildMomentsInvR() {
    const auto& [obsTri, srcTri] = getTriPair();

    momentsInvR = { 0.0, vec3cd::Zero(), 0.0 };
    for (const auto& [obs, obsWeight] : obsTri.quads) {
        const auto& [scaRad, vecRad] = srcTri.getIntegratedInvR(obs);

        auto& [m0, m1, m11] = momentsInvR;
        m0 += scaRad * obsWeight;
        m1 += vecRad * obsWeight;
        m11 += obs.dot(vecRad) * obsWeight;
    }

    // TODO: Symmetric case
}
*/