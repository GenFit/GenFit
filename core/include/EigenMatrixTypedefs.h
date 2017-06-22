#pragma once

#include <Eigen/Dense>

namespace genfit {
    typedef double Precision;

    typedef Precision Scalar;

    typedef Eigen::Matrix<Precision, 3, 1> Vector3;
    typedef Eigen::Matrix<Precision, 5, 1> Vector5;
    typedef Eigen::Matrix<Precision, 6, 1> Vector6;
    typedef Eigen::Matrix<Precision, 7, 1> Vector7;

    typedef Eigen::Matrix<Precision, 3, 3> Matrix3x3;
    typedef Eigen::Matrix<Precision, 4, 4> Matrix4x4;
    typedef Eigen::Matrix<Precision, 5, 5> Matrix5x5;
    typedef Eigen::Matrix<Precision, 6, 6> Matrix6x6;
    typedef Eigen::Matrix<Precision, 7, 7> Matrix7x7;

    typedef Eigen::Matrix<Precision, Eigen::Dynamic, 1> VectorDynamic;
    typedef Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic> MatrixDynamic;

    typedef Matrix3x3 Matrix3x3Sym;
    typedef Matrix4x4 Matrix4x4Sym;
    typedef Matrix5x5 Matrix5x5Sym;
    typedef Matrix6x6 Matrix6x6Sym;
    typedef Matrix7x7 Matrix7x7Sym;
}