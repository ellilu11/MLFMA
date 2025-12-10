#include "rwg.h"

using namespace std;

RWGVec importRWG(
    const filesystem::path& fpath, 
    const std::vector<vec3d>& vertices,
    const TriVec& triangles,
    const shared_ptr<Src> Einc)
{
    ifstream file(fpath);
    string line;
    if (!file) throw std::runtime_error("Unable to find file");
    RWGVec rwgs;

    while (getline(file, line)) {
        istringstream iss(line);
        Eigen::Vector4i idx;

        if (iss >> idx)
            rwgs.emplace_back(make_shared<RWG>(idx, vertices, triangles, Einc));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return rwgs;
}

RWG::RWG(const Eigen::Vector4i& idx,
    const std::vector<vec3d>& vertices,
    const TriVec& triangles,
    const std::shared_ptr<Src> Einc)
    : v0(vertices[idx[0]]), v1(vertices[idx[1]]),
      triPlus(triangles[idx[2]]), triMinus(triangles[idx[3]]),
      vCenter((v0+v1)/2.0),
      Einc(Einc)
{
    for (const auto& vIdx : triPlus->getVidx())
        if (vIdx != idx[0] && vIdx != idx[1])
            vPlus = vertices[vIdx];

    for (const auto& vIdx : triMinus->getVidx())
        if (vIdx != idx[0] && vIdx != idx[1])
            vMinus = vertices[vIdx];

    leng = (v0-v1).norm();

    // buildRHS();

    buildCurrent();

    // std::cout << '(' << v0 << ") (" << v1 << ") ("
    //    << vPlus << ") (" << vMinus << ") " << leng << '\n';
};

void RWG::buildRHS() {

    cmplx rhs(0,0);

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        rhs += weightPlus * exp(iu*kvec.dot(quadNode)) 
                    * (quadNode - vPlus).dot(Einc->pol);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        rhs += weightMinus * exp(iu*kvec.dot(quadNode))
                    * (vMinus - quadNode).dot(Einc->pol);

    rhs *= -Einc->amplitude * leng;
    
}

void RWG::buildCurrent() {
    current = 1.0;

    // TODO: Predict current from rhs vector

}

/* getRadAlongDir(X,kvec)
 * Return the integrated radiated amplitude due to this
 * RWG at cartesian point X along direction kvec
 * X    : observation point
 * kvec : wavevector 
 */ 
vec3cd RWG::getRadAlongDir(
    const vec3d& X, const vec3d& kvec) const {
    
    vec3cd coeff = vec3cd::Zero();

    auto [nodesPlus, weightPlus] = triPlus->getQuads();
    for (const auto& quadNode : nodesPlus)
        coeff += weightPlus * exp(iu*kvec.dot(X-quadNode))
                    * (vPlus - quadNode);

    auto [nodesMinus, weightMinus] = triMinus->getQuads();
    for (const auto& quadNode : nodesMinus)
        coeff += weightMinus * exp(iu*kvec.dot(X-quadNode))
                    * (quadNode - vMinus);

    return current * leng * coeff;
}
