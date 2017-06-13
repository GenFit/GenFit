#pragma once

#include "EigenMatrixTypedefs.h"
#include "RKTools.h"

namespace genfit {

    template <unsigned int rows, unsigned int cols>
    Eigen::Matrix<double, rows, cols> RKMatrixToEigenMatrix(const RKMatrix<rows, cols>& rkMatrix) {
        Eigen::Matrix<double, rows, cols> eigenMatrix;

        for (unsigned int row=0; row < rows; ++row) {
            for (unsigned int col=0; col < cols; ++col) {
                eigenMatrix(row, col) = rkMatrix[cols*row + col];
            }
        }

        return eigenMatrix;
    }


    template <unsigned int rows, unsigned int cols>
    RKMatrix<rows, cols> eigenMatrixToRKMatrix(const Eigen::Matrix<double, rows, cols>& eigenMatrix) {
        RKMatrix<rows, cols> rkMatrix;

        for (unsigned int row=0; row < rows; ++row) {
            for (unsigned int col=0; col < cols; ++col) {
                rkMatrix(row, col) = eigenMatrix(row, col);
            }
        }

        return rkMatrix;
    }

}