#include <filesystem>
#include <random>
#include "sources/dipole.h"
#include "sources/rwg.h"

using namespace std; // TODO: Remove

template <class dist0, class dist1 = dist0, class dist2 = dist0>
SrcVec makeDipoles(const Config& config, const shared_ptr<PlaneWave> Einc)
{

    SrcVec dipoles;

    random_device rd;
    mt19937 gen(rd());

    dist0 rand0(0, 1);
    dist1 rand1(0, 1);
    dist2 rand2(0, 1);

    dist0 prand0(0, 1);
    dist1 prand1(0, 1);

    if (config.pdist == Dist::UNIFORM) {
        double lim = config.rootLeng/2.0;
        rand0 = dist0(-lim, lim);
        rand1 = dist1(-lim, lim);
        rand2 = dist2(-lim, lim);
    }

    for (int n = 0; n < config.nsrcs; ++n) {
        vec3d X = [&] {
            double r, th, ph, z;

            switch (config.pdist) {
                case Dist::UNIFORM:
                    return vec3d(rand0(gen), rand1(gen), rand2(gen));

                // 2D Gaussian at z = 0; TODO: 3D Gaussian
                case Dist::GAUSSIAN: {
                    r = sqrt(-2.0 * log(rand0(gen)));
                    th = PI / 2.0;
                    ph = 2.0 * PI * rand1(gen);
                    return Math::fromSph(vec3d(r, th, ph));
                }

                case Dist::SPHERE: {
                    r = 0.90 * (config.rootLeng / 2.0);
                    th = acos(2.0*rand0(gen) - 1.0);
                    ph = 2.0 * PI * rand1(gen);
                    return Math::fromSph(vec3d(r, th, ph));
                }

                case Dist::CYLINDER: {
                    r = 0.45 * (config.rootLeng / 2.0);
                    ph = 2.0 * PI * rand1(gen);
                    z = 0.90 * config.rootLeng * (rand1(gen) - 0.5);
                    return Math::fromCyl(vec3d(r, ph, z));
                }
            }
            }();

        vec3d P = [&] {
            double th, ph;

            switch (config.qdist) {
                case QDist::UNIFORM:
                    return vec3d(Phys::p0, 0.0, 0.0);

                case QDist::RANDSIGN: {
                    uniform_int_distribution randi(0, 1);
                    return vec3d(Math::sign(randi(gen))*Phys::p0, 0.0, 0.0);
                }

                case QDist::RANDOM: {
                    th = acos(2.0*prand0(gen) - 1.0);
                    ph = 2.0 * PI * prand1(gen);
                    return Math::fromSph(vec3d(Phys::p0, th, ph));
                }
            }
            }();

        dipoles.push_back(make_shared<Dipole>(Einc, X, P));
    }

    return dipoles;
}

SrcVec importDipoles(
    const filesystem::path& fpath,
    const shared_ptr<PlaneWave>& Einc)
{
    ifstream inFile(fpath);
    if (!inFile) throw runtime_error("Unable to find file");
    string line;
    SrcVec dipoles;

    while (getline(inFile, line)) {
        istringstream iss(line);

        vec3d pos, dip;
        if (iss >> pos >> dip)
            dipoles.emplace_back(make_shared<Dipole>(Einc, pos, dip));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return dipoles;
}

vector<vec3d> importVertices(const filesystem::path& fpath) {
    ifstream file(fpath);
    if (!file) throw std::runtime_error("Unable to find file");
    string line;
    vector<vec3d> vList;

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
    const filesystem::path& fpath, const vector<vec3d>& vList, const Precision prec)
{
    ifstream file(fpath);
    string line;
    if (!file) throw std::runtime_error("Unable to find file");
    TriVec triangles;

    while (getline(file, line)) {
        istringstream iss(line);
        vec3i vIdx;

        if (iss >> vIdx)
            triangles.emplace_back(make_shared<Triangle>(vIdx, vList, prec));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return triangles;
}

SrcVec importRWG(
    const filesystem::path& vpath,
    const filesystem::path& tpath,
    const filesystem::path& rpath,
    const Precision quadPrec,
    const shared_ptr<PlaneWave> Einc)
{
    auto vertices = importVertices(vpath);

    auto triangles = importTriangles(tpath, vertices, quadPrec);

    ifstream file(rpath);
    string line;
    if (!file) throw std::runtime_error("Unable to find file");
    SrcVec rwgs;

    while (getline(file, line)) {
        istringstream iss(line);
        Eigen::Vector4i idx;

        if (iss >> idx)
            rwgs.emplace_back(make_shared<RWG>(Einc, idx, vertices, triangles));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return rwgs;
}

std::filesystem::path makePath(const Config& config) {
    std::string distStr =
        [&]() -> std::string {
        switch (config.pdist) {
            case Dist::UNIFORM:    return "uniform";
            case Dist::GAUSSIAN:   return "gauss";
            case Dist::SPHERE:     return "sphere";
            case Dist::CYLINDER:   return "cyl";
        }
        }();

    return
        std::filesystem::path("config") /
        (distStr + "_n" + std::to_string(config.nsrcs) + ".txt");
}

pair<SrcVec, shared_ptr<PlaneWave>> importFromConfig(const Config& config) 
{
    cout << " Importing config...\n";

    const auto fpath = makePath(config);

    shared_ptr<PlaneWave> Einc = make_shared<PlaneWave>();

    /* Dipole sources
    SrcVec srcs;
    switch (config.mode) {
        case Mode::READ:
            srcs = importDipoles(fpath, Einc);
            break;
            
        case Mode::WRITE: {
            srcs = makeDipoles<uniform_real_distribution<double>>(config, Einc);

            ofstream srcFile(fpath);
            for (const auto& src : srcs) srcFile << *(dynamic_pointer_cast<Dipole>(src));
            break;
        }
    }
    cout << "   Source file:     " << fpath.generic_string() << '\n';
    */

    // RWG sources
    const string configPath = "config/n"+to_string(config.nsrcs)+"/";
    auto srcs = importRWG(configPath+"vertices.txt",
                          configPath+"faces.txt",
                          configPath+"rwgs.txt",
                          config.quadPrec,
                          Einc);
    //
    
    cout << fixed << setprecision(3);
    cout << "   Mode:            " << (config.mode == Mode::READ ? "READ" : "WRITE") << '\n';
    cout << "   # Sources:       " << srcs.size() << '\n';
    // cout << "   RWG quad rule:   " << Triangle::prec2Int(config.quadPrec) << "-point\n";
    cout << "   Digit precision: " << config.digits << '\n';
    cout << "   Interp order:    " << config.interpOrder << '\n';
    cout << "   Max node RWGs:   " << config.maxNodeSrcs << '\n';
    cout << "   Root length:     " << config.rootLeng << '\n';
    cout << "   Wave number:     " << Einc->wavenum << "\n\n";

    return make_pair(srcs, Einc);
}

void printSols(SrcVec& srcs, const std::string& fname)
{
    namespace fs = std::filesystem;
    fs::path dir = "out/sol";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    ofstream file(dir/fname);

    file << setprecision(15) << scientific;

    for (const auto& src : srcs)
        src->printSol(file);

}