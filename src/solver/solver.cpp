#include "solver.h"
#include "../mesh/rwg.h"
#include "../mesh/triangle.h"

vecXcd Solver::currents;
vecXcd Solver::lvec;
vecXcd Solver::rvec;

Solver::Solver(const SrcVec& srcs, std::unique_ptr<FMM::Nearfield> nf)
    : nf(std::move(nf)), nsols(srcs.size())
{
    currents = vecXcd::Zero(nsols);
    lvec = vecXcd::Zero(nsols);
    rvec = vecXcd::Zero(nsols);

    // lvec = r = ZI - w = -w assuming I = 0 initially
    for (const auto& src : srcs)
        lvec[src->getIdx()] = -src->getVoltage();
};

void Solver::printSols(const std::string& fname, const vecXcd& sols) {
    std::filesystem::path dir = "out/sol";
    makeDir(dir);
    std::ofstream file(dir/fname);

    file << std::setprecision(15) << std::scientific;
    for (const auto& sol : sols) file << sol << '\n';
}

/* printScattered()
 * Compute and print the scattered far field in the direction of a great circle
 * centered at the origin, with angles theta in [0, pi] and phi = 0 or pi
 * (fold back to [0, pi] to get full coverage of great circle)
 */
void Solver::printScattered(const SrcVec& srcs,
    const std::filesystem::path& dir, const std::string& fname, int nangles)
{
    makeDir(dir);
    std::ofstream farfile(dir/fname), thfile(dir/"angles.txt");
    farfile << std::setprecision(15) << std::scientific;

    std::cout << " Computing scattered farfield...\n";

    double E0 = Exct::Eincs[0]->amplitude;
    double k = config.k, k2 = k*k;
    double rcsSum = 0.0;
    for (int ith = 0; ith < nangles; ++ith) { // 2*nth to cover great circle
        double alpha = (ith+0.5)*PI/static_cast<double>(nangles); // in [0, 2*pi]
        double theta = (alpha < PI) ? alpha : 2*PI - alpha; // fold back to [0, pi]
        double phi = (alpha < PI) ? 0.0 : PI; // fold back to [0, 2pi]
        assert(theta >= 0.0 && theta <= PI);
    //for (int iph = 0; iph < nangles; ++iph) {
    //    double theta = PI / 2.0;
    //    double phi = (iph+0.5)*PI/static_cast<double>(nangles); // in [0, pi]

        vec3d rhat = Math::fromSph(vec3d(1.0, theta, phi));
        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += currents[src->getIdx()] * src->getFarAlongDir(k*rhat);

        // Get theta and phi components of scattered far field
        vec2cd far = iu*k*Phys::eta * Math::toThPh(theta, phi) * dirFar;
        double rcs = 4.0*PI * far.squaredNorm() / (E0*E0);

        rcsSum += rcs;
        farfile << rcs << '\n'; 
        thfile << alpha << ' ' << theta << ' ' << phi << '\n';
        // thfile << phi << ' ' << phi << ' ' << theta << '\n';
    }

    std::cout << "   Mean RCS: "
        << std::setprecision(9) << rcsSum/nangles << std::setprecision(3) << "\n";
    std::cout << "   File: " << (dir/fname).generic_string() << "\n\n";
}

/* printSurfCurrents()
 * Get the surface current at every triangle by summing contributions from adjacent RWGs
 * Use the center of the triangle as the evaluation point for the RWG function
 */
void Solver::printSurfCurrents(const SrcVec& srcs,
    const std::filesystem::path& dir, const std::string& fname) 
{
    makeDir(dir);
    std::ofstream outfile(dir/fname);

    for (const auto& tri : Mesh::glTris) {
        vec3cd J = vec3cd::Zero();
        vec3d center = tri.getCenter();
        const auto& triToRWG = Mesh::glTriToRWGs[tri.getIdx()];

        for (int i = 0; i < 3; ++i) {
            auto [iSrc, isMinus] = triToRWG[i];
            auto rwg = dynamic_pointer_cast<Mesh::RWG>(srcs[iSrc]);

            J += currents[iSrc] * rwg->evaluate(center, isMinus);
        }

        std::cout << center[0] << ' ' << center[1] << ' ' << center[2] << ' ' 
            << J.norm() << '\n';
    }
}