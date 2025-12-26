#pragma once

#include "source.h"
#include "triangle.h"

class RWG final : public Source {
public:

    RWG(std::shared_ptr<PlaneWave>,
        const Eigen::Vector4i&, 
        const std::vector<vec3d>&,
        const TriVec&);

    double getLeng() const { return leng; }

    void buildVoltage() override;

    void buildCurrent() override;

    vec3d getCenter() const override { return center; } 

    vec3cd getRadAlongDir(const vec3d&, const vec3d&, bool = 0) const override;

    // vec3cd getRadAtPoint(const vec3d&) const override;

    cmplx getIntegratedRad(const std::shared_ptr<Source>) const override;

private:
    std::array<std::shared_ptr<Triangle>, 2> tris;
    std::array<vec3d, 2> Xpm; // Non-common vertices

    vec3d X0; // 1st common vertex
    vec3d X1; // 2nd common vertex

    vec3d center; // midpoint of common edge
    double leng; // length of common edge

};