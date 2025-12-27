#include "solver.h"

Solver::Solver(const SrcVec& srcs)
    : nsols(srcs.size()),
      qvec(vecXcd::Zero(nsols)),
      sols(vecXcd::Zero(nsols))

{
    // Guess zero initial current for all sources
    for (int idx = 0; idx < nsols; ++idx)
        qvec[idx] = -srcs[idx]->getVoltage();

    qvec.normalize();

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