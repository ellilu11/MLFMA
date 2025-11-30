#include "triangle.h"

using namespace std;

std::vector<vec3d> importVertices(const filesystem::path& fpath) {
    ifstream file(fpath);
    if (!file) throw std::runtime_error("Unable to find file");
    string line;
    std::vector<vec3d> vList;

    while (getline(file, line)) {
        istringstream iss(line);
        vec3d vertex;

        if (iss >> vertex)
            vList.push_back(vertex);
        else
            throw std::runtime_error("Unable to parse line");
    }

    return vList;
}

TriVec importTriangles(
    const filesystem::path& fpath, const std::vector<vec3d>& vList) 
{
    ifstream file(fpath);
    string line;
    if (!file) throw std::runtime_error("Unable to find file");
    TriVec triangles;

    while (getline(file, line)) {
        istringstream iss(line);
        vec3i vIdx;

        if (iss >> vIdx)
            triangles.emplace_back(make_shared<Triangle>(vIdx,vList));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return triangles;
}

void Triangle::buildQuadNodes(const Precision prec) {

    switch (prec) {
        case Precision::LOW :
            quadNodes.push_back(center);
            quadWeight = 1.0/2.0;
            break;

        case Precision::MEDIUM :
            quadNodes.resize(3);
            quadNodes[0] = 2.0/3.0*vertices[0] + 1.0/6.0*vertices[1] + 1.0/6.0*vertices[2];
            quadNodes[1] = 1.0/6.0*vertices[0] + 2.0/3.0*vertices[1] + 1.0/6.0*vertices[2];
            quadNodes[2] = 1.0/6.0*vertices[0] + 1.0/6.0*vertices[1] + 2.0/3.0*vertices[2];
            quadWeight = 1.0/6.0;
            break;

        case Precision::HIGH :
            quadNodes.resize(7);
            // TODO: 7-point quadrature
            break;
    }
}