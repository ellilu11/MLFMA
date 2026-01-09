#pragma once

#include <Eigen/Dense>

using cmplx = std::complex<double>;
using realVec = std::vector<double>;
using cmplxVec = std::vector<cmplx>;

using pair2i = std::pair<int, int>;
using pair2d = std::pair<double, double>;
using pair2cd = std::pair<cmplx, cmplx>;

using vec3i = Eigen::Vector3i;
using vec3d = Eigen::Vector3d;
using vecXd = Eigen::VectorXd;

using vec2cd = Eigen::Vector2cd;
using vec3cd = Eigen::Vector3cd;
using vecXcd = Eigen::VectorXcd;

using mat3d = Eigen::Matrix3d;
using mat23d = Eigen::Matrix<double, 2, 3>;
using matXcd = Eigen::MatrixXcd;

using interpPair = std::pair<vecXd, int>;
using quadPair = std::pair<vec3d, double>;

template <typename T>
std::vector<T> operator+ (const std::vector<T>& zs, const std::vector<T>& ws) {
    assert(zs.size() == ws.size());

    std::vector<T> sum;
    for (size_t i = 0; i < zs.size(); ++i)
        sum.push_back(zs[i] + ws[i]);
    return sum;
}

/*std::ostream& operator<< (std::ostream& os, cmplx z) {
    //char sign = z.imag() >= 0.0 ? '+' : '-';
    //os << z.real() << sign << abs(z.imag()) << 'i';

    os << z.real() << ' ' << z.imag();

    return os;
}*/

std::ostream& operator<< (std::ostream& os, const vec3d& X) {
    os << X[0] << " " << X[1] << " " << X[2];
    return os;
}

std::ostream& operator<< (std::ostream& os, const vec2cd& X) {
    os << X[0].real() << " " << X[0].imag() << " " 
       << X[1].real() << " " << X[1].imag();
    return os;
}

std::ostream& operator<< (std::ostream& os, const vec3cd& X) {
    os << X[0].real() << " " << X[0].imag() << " "
       << X[1].real() << " " << X[1].imag() << " "
       << X[2].real() << " " << X[2].imag();
    return os;
}

std::ostream& operator<< (std::ostream& os, const Eigen::Array3d& X) {
    os << X[0] << " " << X[1] << " " << X[2];
    return os;
}

std::ostream& operator<< (std::ostream& os, const vecXcd& X) {
    for (int i = 0; i < X.rows(); ++i)
        os << X[i] << "\n";

    return os;
}

std::istream& operator>>(std::istream& is, vec3d& X) {
    double x, y, z;
    if (is >> x >> y >> z)
        X = vec3d{ x, y, z };
    return is;
}

std::istream& operator>>(std::istream& is, vec3i& X) {
    int x, y, z;
    if (is >> x >> y >> z)
        X = vec3i{ x, y, z };
    return is;
}

std::istream& operator>>(std::istream& is, Eigen::Vector4i& X) {
    int x, y, z, w;
    if (is >> x >> y >> z >> w)
        X = Eigen::Vector4i{ x, y, z, w };
    return is;
}

/*
template <typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& X) {
    T x, y, z;
    if (is >> x >> y >> z)
        X = std::vector<T>{ x, y, z };
    return is;
}*/

std::array<bool,3> operator> (const vec3d& x, const vec3d& y) {
    std::array<bool,3> bools{ x[0] > y[0], x[1] > y[1], x[2] > y[2] };
    return bools;
}