#pragma once

#include <filesystem>
#include <random>
#include "math.h"

class Triangle;

using TriVec = std::vector<std::shared_ptr<Triangle>>;

class Triangle {
public:
    Triangle()
        : vIdx(vec3i::Zero()),
          vertices( {zeroVec,zeroVec,zeroVec}), 
          center(zeroVec), area(0.0) { };

    Triangle(const vec3i& vIdx, const std::vector<vec3d>& vList) 
        : vIdx(vIdx),
          vertices( { vList[vIdx[0]], vList[vIdx[1]], vList[vIdx[2]] } )
          // center( (vertices[0] + vertices[1] + vertices[2]) / 3.0 )
    {
    };

    vec3i getVidx() { return vIdx; }

    std::array<vec3d, 3> getVertices() { return vertices; }

    vec3d getCenter() const { return center; }

    double getArea() const { return area; }

private:
    vec3i vIdx;
    std::array<vec3d,3> vertices;
    vec3d center;
    double area;
};