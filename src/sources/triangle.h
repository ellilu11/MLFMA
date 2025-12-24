#pragma once

class Triangle;

using TriVec = std::vector<std::shared_ptr<Triangle>>;

class Triangle {
public:
    friend class RWG;

    Triangle()
        : vIdx(vec3i::Zero()),
          Xs( {zeroVec,zeroVec,zeroVec}) 
          {};

    Triangle(
        const vec3i&,
        const std::vector<vec3d>&,
        const Precision);

    // vec3i getVidx() { return vIdx; }

    //std::array<vec3d, 3> getVertices() { return Xs; }

    //std::array<vec3d, 3> getDistVecs() { return Ds; }

    std::pair<std::vector<vec3d>, double> getQuads() {
        return std::make_pair(quadNodes, quadWeight);
    }

    static int prec2Int(const Precision);

    void buildQuads(const Precision);

private:
    std::array<vec3d,3> Xs; // vertices
    std::array<vec3d,3> Dps; // Dps[i] = Xs[i+1] - Xs[i]
    std::array<vec3d,3> Dms; // Dms[i] = Xs[i] - Xs[i-1]

    vec3i vIdx;

    std::vector<vec3d> quadNodes;
    double quadWeight;

    vec3d nhat; // surface normal unit vector

    // double area; // not needed due to cancellation
};