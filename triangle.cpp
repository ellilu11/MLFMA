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