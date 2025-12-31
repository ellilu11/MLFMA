#include "solver.h"

Solver::Solver(const SrcVec& srcs, std::shared_ptr<Node> root, int numIter, double EPS)
    : root(std::move(root)),
      numSols(srcs.size()),
      numIter(numIter),
      EPS(EPS),
      Qmat(matXcd(numSols, 1)),
      gvec(vecXcd::Zero(numIter+1)),
      currents(vecXcd::Zero(numSols)), // Guess zero currents
      qvec(std::make_shared<vecXcd>(vecXcd::Zero(numSols))),
      sols(std::make_shared<vecXcd>(vecXcd::Zero(numSols)))
{
    // Guess zero currents
    // std::transform
    for (int idx = 0; idx < numSols; ++idx)
        (*qvec)[idx] = -srcs[idx]->getVoltage();

    g0 = (*qvec).norm();
    gvec[0] = g0;

    (*qvec).normalize(); // qvec_0
    Qmat.col(0) = *qvec; // store qvec as first column of Qmat
}

void Solver::updateSols(int k) {
    root->buildMpoleCoeffs();

    root->buildLocalCoeffs();

    Leaf::evaluateSols();
}

void Solver::iterateArnoldi(int k) {
    assert(Qmat.cols() == k+1);

    // Do Arnoldi iteration
    vecXcd hcol(k+2);
    for (int i = 0; i <= k; ++i) {
        const auto& q_i = Qmat.col(i); // get qvec_i
        cmplx h = q_i.dot(*sols); // Hermitian dot
        *sols -= h * q_i;
        hcol[i] = h;
    }
    hcol[k+1] = (*sols).norm();

    // Replace present qvec with new qvec
    *qvec = (*sols) / hcol[k+1]; // .normalized();

    // Store new qvec as new column of Qmat
    Qmat.conservativeResize(numSols, k+2);
    Qmat.col(k+1) = *qvec;

    // Store new hcol as new column of Hmat
    Hmat.conservativeResize(k+2, k+1); // 2x1, 3x2, 4x3 ...
    Hmat.row(k+1).setZero();
    Hmat.col(k) = hcol;
}

void Solver::applyGivensRotation(
    vecXcd& hcol, vecXcd& vcos, vecXcd& vsin, int k) {

    assert(hcol.size() == k+2);

    for (int i = 0; i < k; ++i) {
        const cmplx hprev = vcos[i] * hcol[i] + vsin[i] * hcol[i+1] ;
        hcol[i+1] = -vsin[i] * hcol[i] + vcos[i] * hcol[i+1];
        hcol[i] = hprev;
    }

    auto [vcos_k, vsin_k] = Math::givensRotation(hcol[k], hcol[k+1]);

    hcol[k] = vcos_k * hcol[k] + vsin_k * hcol[k+1];
    hcol[k+1] = 0.0;

    vcos[k] = vcos_k;
    vsin[k] = vsin_k;
}

void Solver::updateGvec(vecXcd& vcos, vecXcd& vsin, int k) {

    vecXcd hcol_k = Hmat.col(k);
    applyGivensRotation(hcol_k, vcos, vsin, k);
    Hmat.col(k) = hcol_k;

    gvec[k+1] = -vsin[k] * gvec[k];
    gvec[k] = vcos[k] * gvec[k];
}

void Solver::solve() {

    vecXcd vcos = vecXcd::Zero(numIter);
    vecXcd vsin = vecXcd::Zero(numIter);

    int iter = 0;
    do {
        std::cout << " Do iteration #" << iter << '\n';
        auto iter_start = Clock::now();

        updateSols(iter);

        iterateArnoldi(iter);

        updateGvec(vcos, vsin, iter);
        
        resetSols();

        Time fmm_duration_ms = Clock::now() - iter_start;
        std::cout << "   Elapsed time: " << fmm_duration_ms.count() << " ms\n";

    } while (abs(gvec[++iter])/g0 > EPS);

    std::cout << " Solving for current...\n";

    const auto& Hp = Hmat.block(0, 0, Hmat.rows()-1, Hmat.cols());
    vecXcd yvec = Hp.lu().solve(gvec.segment(0, iter));
    currents = Qmat.leftCols(iter) * yvec;

    printSols("solDir.txt");
}

/*void Solver::solve() {

    vecXcd vcos = vecXcd::Zero(numIter);
    vecXcd vsin = vecXcd::Zero(numIter);

    // std::string outStr = "out/sol/sol.txt";
    // std::filesystem::remove(outStr);
    // std::ofstream outFile(outStr);

    int iter = 0;
    for (int iter = 0; iter < numIter; ++iter) {
        std::cout << " Do iteration #" << iter << '\n';
        auto iter_start = Clock::now();

        updateSols(iter);

        if (iter == numIter-1) printSols("sol.txt");

        resetSols();

        Time fmm_duration_ms = Clock::now() - iter_start;
        std::cout << "   Elapsed time: " << fmm_duration_ms.count() << " ms\n";

    }
}*/

/*void Solver::solve() {

    //
    namespace fs = std::filesystem;
    fs::path dir = "out/sol";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::filesystem::remove(dir/"sol.txt");
    std::ofstream file(dir/"sol.txt", std::ios::app);
    //

    vecXcd vcos = vecXcd::Zero(numIter);
    vecXcd vsin = vecXcd::Zero(numIter);

    for (int iter = 0; iter < numIter; ++iter) {
        updateSols(iter);

        iterateArnoldi(iter);

        updateGvec(vcos, vsin, iter);

        resetSols();

        //
        const auto& Hp = Hmat.block(0, 0, Hmat.rows()-1, Hmat.cols());
        vecXcd yvec = Hp.lu().solve(gvec.segment(0, iter+1));

        currents = Qmat.leftCols(iter+1) * yvec;

        for (const auto& curr : currents)
            os << curr.real() << ' ';
        os << '\n';
        //

        const double err = abs(gvec[iter+1])/g0;
        if (err < EPS) break;
    }
}*/

void Solver::printSols(const std::string& fname) {
    namespace fs = std::filesystem;
    fs::path dir = "out/sol";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream file(dir/fname);

    file << std::setprecision(15) << std::scientific;

    //for (const auto& sol : *sols)
    //    file << sol.real() << ' ' << sol.imag() << '\n';

    for (const auto& curr : currents)
        file << curr.real() << ' ' << curr.imag() << '\n';
}