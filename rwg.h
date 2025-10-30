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

    vec3d getCenter() { return vCenter; }

    cmplx buildExcitation();

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
    double b;
};