#include "angles.h"

void FMM::Angles::buildAngularSamples(int level)
{
    const double wavenum = Node::wavenum;
    const double nodeLeng = Node::config.rootLeng / pow(2.0, level);

    // Use excess bandwidth formula
    const int tau = ceil((1.73*wavenum*nodeLeng +
        2.16*pow(Node::config.digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));

    L = floor(0.50*tau); // TODO: Find optimal formula

    // Construct thetas
    const int nth = tau+1;
    std::tie(thetas, weights) = Interp::gaussLegendre(nth, EPS_NR, 0.0, PI);

    // Absorb sin(theta) into weights
    std::transform(weights.begin(), weights.end(), thetas.begin(), weights.begin(),
        [](double weight, double theta) { return weight * sin(theta); }
    );

    // Construct phis
    const int nph = 2*nth;
    phis.resize(nph);

    for (int iph = 0; iph < nph; ++iph)
        phis[iph] = 2.0*PI*iph/static_cast<double>(nph);

    std::cout << "   (" << level << "," << thetas.size() << "," << phis.size() << ")\n";
}

void FMM::Angles::buildAngularMatrices() {
    const auto [nth, nph] = getNumAngles();
    const size_t nDirs = nth*nph;

    khat.resize(nDirs);
    toThPh.resize(nDirs);
    ImRR.resize(nDirs);

    size_t iDir = 0;
    for (int ith = 0; ith < nth; ++ith) {
        const double theta = thetas[ith];

        for (int iph = 0; iph < nph; ++iph) {
            const double phi = phis[iph];

            khat[iDir] = Math::fromSph(vec3d(1.0, theta, phi));
            toThPh[iDir] = Math::toThPh(theta, phi);
            ImRR[iDir] = Math::ImRR(khat[iDir]);

            ++iDir;
        }
    }
}

void FMM::Angles::printAngles(std::ofstream& thfile, std::ofstream& phfile) {
    thfile << std::setprecision(15);
    phfile << std::setprecision(15);

    for (const auto& theta : thetas)
        thfile << theta << '\n';

    for (const auto& phi : phis)
        phfile << phi << '\n';
}