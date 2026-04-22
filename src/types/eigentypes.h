#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using vec2i = Eigen::Vector2i;
using vec3i = Eigen::Vector3i;
using vec4i = Eigen::Vector4i;
using vec3d = Eigen::Vector3d;
using vecXd = Eigen::VectorXd;
using vec2cd = Eigen::Vector2cd;
using vec3cd = Eigen::Vector3cd;
using vecXcd = Eigen::VectorXcd;

using arrXcd = Eigen::ArrayXcd;

using mat3d = Eigen::Matrix3d;
using mat23d = Eigen::Matrix<double, 2, 3>;
using matXcd = Eigen::MatrixXcd;

template <typename T>
using sparseMat = Eigen::SparseMatrix<T, Eigen::RowMajor>;

vec3i operator> (const vec3d& X, const vec3d& Y) {
    vec3i bools{ X[0] > Y[0], X[1] > Y[1], X[2] > Y[2] };
    return bools;
}

std::ostream& operator<< (std::ostream& os, const vec3d& X) {
    // os << X[0] << " " << X[1] << " " << X[2];
    os << '(' << X[0] << " " << X[1] << " " << X[2] << ')';
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

std::istream& operator>>(std::istream& is, vec4i& X) {
    int x, y, z, w;
    if (is >> x >> y >> z >> w)
        X = vec4i{ x, y, z, w };
    return is;
}