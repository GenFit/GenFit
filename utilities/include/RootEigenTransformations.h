#pragma once

#include <eigen3/Eigen/Dense>


namespace genfit {

    template <unsigned int dim>
    Eigen::Matrix<double, dim, 1> rootVectorToEigenVector(const TVectorD& rootVector) {
        const unsigned int rootVectorRows = rootVector.GetNrows();

        double eigenArray[dim];
        const auto* rootArray = rootVector.GetMatrixArray();

        std::copy(rootArray,
                  rootArray + rootVectorRows,
                  std::begin(eigenArray));

        return Eigen::Map<Eigen::Matrix<double, dim, 1>>(eigenArray);
    }

    template <unsigned int dim>
    TVectorD eigenVectorToRootVector(const Eigen::Matrix<double, dim, 1>& eigenVector) {
        const unsigned int eigenVectorRows = eigenVector.rows();

        const double* eigenArray = eigenVector.data();

        TVectorD rootVector(eigenVectorRows);
        std::copy(eigenArray,
                  eigenArray + eigenVectorRows,
                  rootVector.GetMatrixArray());

        return rootVector;
    }

    template <unsigned int dim>
    Eigen::Matrix<double, dim, dim> rootMatrixSymToEigenMatrix(const TMatrixDSym& rootMatrix) {
        Eigen::Matrix<double, dim, dim> eigenMatrix;

        for (unsigned int row=0; row<dim; ++row) {
            for (unsigned int col=0; col<dim; ++col) {
                eigenMatrix(row, col) = rootMatrix(row, col);
                eigenMatrix(col, row) = eigenMatrix(row, col);
            }
        }

        return eigenMatrix;
    }

    template <unsigned int dim>
    TMatrixDSym eigenMatrixToRootMatrixSym(const Eigen::Matrix<double, dim, dim>& eigenMatrix) {
        TMatrixDSym rootMatrix(dim);

        for (unsigned int row=0; row<dim; ++row) {
            for (unsigned int col=0; col<dim; ++col) {
                rootMatrix(row, col) = eigenMatrix(row, col);
                rootMatrix(col, row) = rootMatrix(row, col);
            }
        }

        return rootMatrix;
    }

}
