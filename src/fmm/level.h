#pragma once

#include "../maps.h"
#include "../math.h"
#include "../mesh/mesh.h"
#include "fmm.h"

struct FMM::Level {

public:
    Level() = default;

    Level(int lvl) : lvl(lvl)
    {
        buildAngularSamples();
        buildAngularMatrices();
        buildTranslationTable();

        /*std::cout << "   (" << lvl << ","
            << thetas.size() << "," << phis.size() << ") ";*/
    }

    void buildInterpTables(const Level& srcLevel) {
        assert(lvl < maxLevel);
        buildInterpTheta(srcLevel);
        buildInterpPhi(srcLevel);
    }

    static void buildDists();

    static void clearDists() {
        dists = std::vector<double>();
        rhats = std::vector<vec3d>();
        dXs = std::vector<vec3d>();
    }

    pair2i getNumAngles() const {
        return std::make_pair(thetas.size(), phis.size());
    }

    size_t getNumDirs() const {
        return thetas.size() * phis.size();
    }

    void printAngles(std::ofstream& thfile, std::ofstream& phfile) {
        thfile << std::setprecision(15);
        for (const auto& theta : thetas) thfile << theta << '\n';

        phfile << std::setprecision(15);
        for (const auto& phi : phis) phfile << phi << '\n';
    }

private:
    void buildAngularSamples();

    void buildAngularMatrices();

    void buildInterpTheta(const Level&);

    void buildInterpPhi(const Level&);

    Map<std::vector<cmplx>> getAlpha();

    HashMap<interpPair> getInterpPsi();

    void buildTranslationTable();

public:
    // M2M interpolation tables
    std::vector<interpPair> interpTheta;
    std::vector<interpPair> interpPhi;

    // M2L translation table
    VecHashMap<arrXcd> transl;

    // Angular quantities
    std::vector<vec3d> khat;
    std::vector<mat23d> toThPh;
    std::vector<mat3d> ImRR;

    std::vector<double> thetas;  // theta samples
    std::vector<double> weights; // weights of theta samples
    std::vector<double> phis;    // phi samples

private:
    static std::vector<double> dists; // unique distances between interacting nodes
    static std::vector<vec3d> rhats;  // unique unit distance vectors between interacting nodes
    static std::vector<vec3d> dXs;    // unique distance vectors between interacting nodes

    int L;     // M2L series truncation number
    int lvl;   // level of this FMM level, with root at level 0
};
