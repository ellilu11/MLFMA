#include "solver.h"

Solver::Solver(const SrcVec& srcs)
    : nsols(srcs.size()),
      Q(matXcd(nsols, 1)),
      qvec(vecXcd::Zero(nsols)),
      sols(vecXcd::Zero(nsols))
{
    // Guess zero initial current for all sources
    for (int idx = 0; idx < nsols; ++idx)
        qvec[idx] = -srcs[idx]->getVoltage();

    qvec.normalize(); // qvec_0
    Q.col(0) = qvec; // store qvec as first column of Q

    // gvec.push_back(qvec.norm());
}

// Assume iter starts from 1 for easier bookkeeping
void Solver::iterateArnoldi(int iter) {

    assert(Q.cols() == iter);

    vecXcd hvec(iter+1);

    // Do Arnoldi iteration
    for (int i = 0; i < iter; ++i) {
        const auto& q_i = Q.col(i); // get qvec_i
        cmplx h = q_i.dot(sols);
        sols -= h * q_i;
        hvec[i] = h;
        // Normalize sols at every Arnoldi iteration?
    }
    hvec[iter] = sols.norm();

    // Replace current qvec with new qvec
    qvec = sols.normalized();

    // Store new qvec as new column of Q
    Q.conservativeResize(nsols, iter+1);
    Q.col(iter) = qvec;

    // Store new hvec as new column of H
    assert(hvec.rows() == iter+1);
    H.conservativeResize(iter+1, iter);
    H.col(iter) = hvec;

    // Reset sols for next iteration
    resetSols();
}

pair2cd Solver::applyGivensRotation(
    vecXcd& hvec, const vecXcd& vcos, const vecXcd& vsin, int k) {

    for (int i = 0; i < k-1; ++i) {
        const cmplx htemp = hvec(i) * vcos(i) + hvec(i+1) * vsin(i) ;
        hvec(i+1) = -hvec(i) * vsin(i) + hvec(i+1) * vcos(i);
        hvec(i) = htemp;
    }

    auto [vcos_k, vsin_k] = Math::givensRotation(hvec(k), hvec(k+1));

    hvec(k) = vcos_k * hvec(k) + vsin_k * hvec(k+1);
    hvec(k+1) = 0.0; // vecXcd::Zero(k+2)

    return make_pair(vcos_k, vsin_k);
}

void Solver::solve(int numIter) {

    const vecXcd vcos = vecXcd::Zero(numIter);
    const vecXcd vsin = vecXcd::Zero(numIter);
    

    for (int iter = 0; iter < numIter; ++iter) {
        // doMLFMA();

        iterateArnoldi(iter);

        vecXcd H_k = H.col(iter);

        auto [vcos_k, vsin_k] = applyGivensRotation(H_k, vcos, vsin, iter);

        // Finish
    }
}

void Solver::printSols(const std::string& fname)
{
    namespace fs = std::filesystem;
    fs::path dir = "out/sol";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream file(dir/fname);

    file << std::setprecision(15) << std::scientific;

    for (const auto& sol : sols)
        file << sol.real() << ' ' << sol.imag() << '\n';

}