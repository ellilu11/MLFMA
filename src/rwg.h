#pragma once

#include <algorithm>
#include <filesystem>
#include <random>
#include "math.h"
#include "src.h"
#include "triangle.h"

class RWG;

using RWGVec = std::vector<std::shared_ptr<RWG>>;

class RWG {
public:
    RWG(const Eigen::Vector4i&, 
        const std::vector<vec3d>&,
        const TriVec&,
        const std::shared_ptr<Src>);

    std::shared_ptr<Triangle> getTriPlus() const { return triPlus; }
    std::shared_ptr<Triangle> getTriMinus() const { return triMinus; }

    vec3d getCenter() const { return vCenter; }

    vec3d getVplus() const { return vPlus; }

    vec3d getVminus() const { return vMinus; }

    double getLeng() const { return leng; }

    cmplx getCurrent() const { return current; }

    cmplx buildRHS();

    cmplx buildCurrent();

private:
    vec3d v0;
    vec3d v1;
    std::shared_ptr<Triangle> triPlus;
    std::shared_ptr<Triangle> triMinus;

    vec3d vCenter;
    vec3d vPlus;
    vec3d vMinus;

    double leng;

    std::shared_ptr<Src> Einc;
    cmplx rhs;
    cmplx current;

};