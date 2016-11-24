#pragma once

#include <eigen3/Eigen/Dense>


namespace genfit {

    Eigen::Matrix<double, 5, 1> rootVectorToEigenVector(const TVectorD& rootVector) {
        const unsigned int rootVectorRows = rootVector.GetNrows();

        double eigenArray[5];
        const auto* rootArray = rootVector.GetMatrixArray();

        std::copy(rootArray,
                  rootArray + rootVectorRows,
                  std::begin(eigenArray));

        return Eigen::Map<Eigen::Matrix<double, 5, 1>>(eigenArray);
    }

    TVectorD eigenVectorToRootVector(const Eigen::Matrix<double, 5, 1>& eigenVector) {
        const unsigned int eigenVectorRows = eigenVector.rows();

        const double* eigenArray = eigenVector.data();

        TVectorD rootVector(eigenVectorRows);
        std::copy(eigenArray,
                  eigenArray + eigenVectorRows,
                  rootVector.GetMatrixArray());

        return rootVector;
    }

    Eigen::Matrix<double, 5, 5> rootMatrixSymToEigenMatrix(const TMatrixDSym& rootMatrix) {
        const unsigned int rootMatrixRows = rootMatrix.GetNrows();
        const unsigned int rootMatrixCols = rootMatrix.GetNcols();

        Eigen::Matrix<double, 5, 5> eigenMatrix;

        for (unsigned int row=0; row<rootMatrixRows; ++row) {
            for (unsigned int col=0; col<rootMatrixCols; ++col) {
                eigenMatrix(row, col) = rootMatrix(row, col);
            }
        }

        return eigenMatrix;
    }

    TMatrixDSym eigenMatrixToRootMatrixSym(const Eigen::Matrix<double, 5, 5>& eigenMatrix) {
        const unsigned int eigenMatrixRows = eigenMatrix.rows();
        const unsigned int eigenMatrixCols = eigenMatrix.cols();

        TMatrixDSym rootMatrix(5);

        for (unsigned int row=0; row<eigenMatrixRows; ++row) {
            for (unsigned int col=0; col<eigenMatrixCols; ++col) {
                rootMatrix(row, col) = eigenMatrix(row, col);
            }
        }

        return rootMatrix;
    }

}
