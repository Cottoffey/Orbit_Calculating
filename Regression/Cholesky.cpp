#include <cmath>
#include "Cholesky.h"


matrix Cholesky::CholeskyDecomp(const matrix &A) {
    matrix L;
    int n = A.size();
    L.resize(n);
    for (int i = 0; i < n; ++i) {
        L[i].resize(n);
        double sum;
        for (int j = 0; j < i; ++j) {
            sum = 0;
            for (int p = 0; p < j; ++p)
                sum += L[i][p] * L[j][p];
            L[i][j] = (A[i][j] - sum) / L[j][j];
        }
        sum = A[i][i];
        for (int p = 0; p < i; ++p)
            sum -= L[i][p] * L[i][p];
        L[i][i] = sqrt(sum);
    }
    return L;
}

vec Cholesky::solveTriangleSLE(const matrix &A, const vec &B, MATRIX_TYPE type) {
    int dim = B.size();
    vec X;
    X.resize(dim);
    if (type == MATRIX_TYPE::LOWER_TRIANGLE) {
        X[0] = B[0] / A[0][0];
        for (int i = 1; i < dim; ++i) {
            double sum = 0;
            for (int k = 0; k < i; ++k)
                sum += A[i][k] * X[k];
            X[i] = (B[i] - sum) / A[i][i];

        }
    } else if (type == MATRIX_TYPE::UPPER_TRIANGLE) {
        X[dim - 1] = B[dim - 1] / A[dim - 1][dim - 1];
        for (int i = dim - 2; i > -1; --i) {
            double sum = 0;
            for (int k = dim - 1; k > i; --k)
                sum += A[i][k] * X[k];
            X[i] = (B[i] - sum) / A[i][i];
        }
    }
    return X;
}

vec Cholesky::CholeskySLE(const matrix &A, const vec &B) {
    matrix L = CholeskyDecomp(A);
    vec y = solveTriangleSLE(L, B, MATRIX_TYPE::LOWER_TRIANGLE);
    return solveTriangleSLE(transpose(L), y, MATRIX_TYPE::UPPER_TRIANGLE);
}

