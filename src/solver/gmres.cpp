#include "gmres.h"
#include "../fmm/farfield.h"
#include "../fmm/nearfield.h"
#include "../fmm/node.h"

GMRES::GMRES(
    const SrcVec& srcs,
    std::unique_ptr<FMM::Nearfield> nf,
    std::unique_ptr<FMM::Farfield> ff,
    std::shared_ptr<FMM::Node> root)
    : Solver(srcs, std::move(nf)),
      ff(std::move(ff)), root(std::move(root)),
      Qmat(matXcd(nsols, 1)),
      gvec(vecXcd::Zero(config.maxIter+1)),
      vcos(vecXcd::Zero(config.maxIter)),
      vsin(vecXcd::Zero(config.maxIter))
{
    doILU = !this->root->isLeaf() && config.maxIter;
    // doILU = config.maxIter;
    if (doILU) {
        buildILU();
        lvec = ilu.solve(lvec); // apply M^(-1) to lvec
    }

    g0 = lvec.norm();   // store g0 for use later
    gvec[0] = g0;       // store g0 as first entry of gvec
    lvec.normalize();   
    Qmat.col(0) = lvec; // store normalized lvec as first column of Qmat
}

/* buildILU()
 * Build ILU decomposition of nearfield matrix for use as preconditioner
 */
void GMRES::buildILU() {
    ilu.setDroptol(config.iluTol);
    ilu.setFillfactor(config.iluFactor);

    std::cout << " Building preconditioner...       ";
    auto start = Clock::now();

    ilu.compute(nf->nearMat);
    
    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
}

/* updateRvec()
 * Update rvec = Z * lvec by evaluating nearfield and farfield contributions
 * Then apply preconditioner M^(-1) to rvec
 */
void GMRES::updateRvec(int k) {
    // FMM main loop: upward pass, downward pass, evaluate sols
    if (!root->isLeaf()) {
        // Upward pass
        ff->buildGlMpoleCoeffs();
        ff->buildMpoleCoeffs(root, true);

        // Downward pass
        ff->glTranslateCoeffs();
        ff->buildLocalCoeffs(root);

        // Evaluate far sols
        ff->evaluateSols();
    }
    // Evaluate near sols
    nf->evaluateSols();

    if (doILU) rvec = ilu.solve(rvec); // apply M^(-1) to rvec = Z * lvec
}

/* iterateArnoldi()
 * Perform one iteration of Arnoldi process to get new lvec, 
 * and update Hmat and Qmat
 */
void GMRES::iterateArnoldi(int k) {
    assert(Qmat.cols() == k+1);

    // Do Arnoldi iteration
    vecXcd hcol(k+2);
    for (int i = 0; i <= k; ++i) {
        const auto& q_i = Qmat.col(i); // get lvec_i
        cmplx h = q_i.dot(rvec); // Hermitian dot
        rvec -= h * q_i;
        hcol[i] = h;
    }
    hcol[k+1] = rvec.norm();

    // Replace present lvec with new lvec
    lvec = rvec / hcol[k+1]; // .normalized();

    // Store new lvec as new column of Qmat
    Qmat.conservativeResize(nsols, k+2);
    Qmat.col(k+1) = lvec;

    // Store new hcol as new column of Hmat
    Hmat.conservativeResize(k+2, k+1); // 2x1, 3x2, 4x3 ...
    Hmat.row(k+1).setZero();
    Hmat.col(k) = hcol;
}

/* applyGivensRotation()
 * Apply Givens rotations to hcol to introduce zeros below diagonal,
 * and store new cos and sin values for next iteration
 */
void GMRES::applyGivensRotation(vecXcd& hcol, int k) {
    assert(hcol.size() == k+2);

    for (int i = 0; i <= k-1; ++i) {
        cmplx hprev = vcos[i] * hcol[i] + vsin[i] * hcol[i+1] ;
        hcol[i+1] = -vsin[i] * hcol[i] + vcos[i] * hcol[i+1];
        hcol[i] = hprev;
    }

    auto [vcos_k, vsin_k] = Math::givensRotation(hcol[k], hcol[k+1]);

    hcol[k] = vcos_k * hcol[k] + vsin_k * hcol[k+1];
    hcol[k+1] = 0.0;

    vcos[k] = vcos_k;
    vsin[k] = vsin_k;
}

/* updateGvec()
 * Update gvec by applying the same Givens rotations to introduce zeros below diagonal,
 * and update gvec values for next iteration
 */
void GMRES::updateGvec(int k) {
    vecXcd hcol_k = Hmat.col(k);
    applyGivensRotation(hcol_k, k);
    Hmat.col(k) = hcol_k;

    gvec[k+1] = -vsin[k] * gvec[k];
    gvec[k] = vcos[k] * gvec[k];
}

/* solve()
 * Main GMRES loop to solve for currents, and print solutions to file
 */
void GMRES::solve(const std::string& fname) {
    // If maxIter = 0, just do one MVM and print rvec
    if (!config.maxIter) {
        updateRvec(0);
        printSols(fname, rvec);
        std::cout << '\n';
        return;
    }

    std::ofstream resfile("out/sol/" + std::string(root->isLeaf() ? "resDir" : "res") + ".txt");
    std::string method = root->isLeaf() ? "direct... " : "FMM...    ";
    std::cout << " Solving for current w/ " << method;
    auto start = Clock::now();

    int iter = 0;
    do {
        if (iter && !(iter%100)) std::cout << " #" << iter << ' ';
        updateRvec(iter);
        iterateArnoldi(iter);
        updateGvec(iter);
        rvec = vecXcd::Zero(nsols); // reset rvec for next iteration
        resfile << abs(gvec[iter+1])/g0 << '\n'; // log the relative residual
    } while (abs(gvec[++iter])/g0 > config.epsIter && iter < config.maxIter); // careful

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n";
    t.printTimes();
    t.resetTimes();

    std::cout << "   # Iterations: " << iter << "\n";

    // Solve least squares problem to get yvec, then compute currents = Qmat * yvec
    const matXcd& Hp = Hmat.block(0, 0, Hmat.rows()-1, Hmat.cols());
    vecXcd yvec = Hp.lu().solve(gvec.segment(0, iter));
    currents = Qmat.leftCols(iter) * yvec;
    std::cout << "   Current norm: " 
        << std::setprecision(9) << currents.norm() << std::setprecision(3) << "\n";

    printSols(fname, currents);
}