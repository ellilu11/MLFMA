#include <fstream>
#include "main.h"
#include "clock.h"
#include "config.h"
#include "fmm/fmm.h"
#include "solver/direct.h"

extern const Config config("config/config.txt");

void mainLoop(const SrcVec& srcs, bool doFMM, bool doIter = true) {
    size_t nsrcs = srcs.size();

    // ==================== Build nodes ========================= //
    auto start = Clock::now();

    bool isRootLeaf = nsrcs <= config.maxNodeSrcs || !doFMM;
    auto root = std::make_shared<FMM::Node>(srcs, 0, nullptr, isRootLeaf);
    root->postProcess();
   
    // ==================== Build nearfield ===================== //
    auto nf = std::make_unique<FMM::Nearfield>(nsrcs);

    // ==================== Build farfield ====================== //
    auto ff = std::make_unique<FMM::Farfield>(root);

    // ==================== Solve for current =================== //
    std::unique_ptr<Solver> solver;
    if (doFMM || (!doFMM && doIter))
        solver = std::make_unique<GMRES>(srcs, std::move(nf), std::move(ff), root);
    else 
        solver = std::make_unique<Direct>(srcs, std::move(nf));

    solver->solve(doFMM ? "sol.txt" : "solDir.txt");

    Time duration_ms = Clock::now() - start;
    std::cout << " " + std::string(doFMM ? "FMM" : "Direct") + " total elapsed time : "
        << duration_ms.count() << " ms\n\n";

    // ==================== Compute scattered field ============= //
    std::string ieStr = getIEStr(config.ie);
    std::transform(ieStr.begin(), ieStr.end(), ieStr.begin(), ::tolower);
    // Solver::printScattered(srcs, "out/ff/sph/px_k1.0z_r5.0_"+ieStr, (doFMM ? "ff_n" : "ffDir_n")+std::to_string(nsrcs)+".txt", 200);
    // Solver::printScattered(srcs, "out/ff/plate/px_k1.0z_plate", std::string(doFMM ? "ff_g" : "ffDir_g")+config.lengStr+".txt", 100);
    Solver::printScattered(srcs, "out/ff/almond/py_k146.6z_almond_"+ieStr, std::string(doFMM ? "ff" : "ffDir")+"_openmp.txt", 500);
    // Solver::printScattered(srcs, "out/ff/glider/px_k41.9z_glider_"+ieStr, std::string(doFMM ? "ff" : "ffDir")+".txt", 500);
}

int main() {
    omp_set_num_threads(config.numThreads);

    importVec<vec3d>("config/pwave.txt", Exct::Eincs);
    // auto srcs = Mesh::importMesh("config/rwg/sph_r5.0/sph_r5.0_n"+std::to_string(config.nsrcs));
    // auto srcs = Mesh::importMesh("config/rwg/rect/rect_g"+config.lengStr+"_n"+std::to_string(config.nsrcs));
    auto srcs = Mesh::importMesh("config/rwg/almond/almond_n"+std::to_string(config.nsrcs));
    // auto srcs = Mesh::importMesh("config/rwg/glider/glider_n"+std::to_string(config.nsrcs));

    constexpr bool doIter = true;
    switch (config.mode) {
        case Mode::FMM: 
            mainLoop(srcs, true); 
            break;
        case Mode::DIR: 
            mainLoop(srcs, false, doIter); 
            break;
        case Mode::FMMDIR: {
            mainLoop(srcs, true);
            FMM::reset();
            mainLoop(srcs, false, doIter);
        } break;
    }

    return 0;
}