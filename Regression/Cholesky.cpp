#include <cmath>
#include "Cholesky.h"


Matrix Cholesky::CholeskyDecomp(const Matrix &A) {
    int n = A.size();
    Matrix L(n);
    for (int i = 0; i < n; i++) {
        L[i].resize(n);
        for (int j = 0; j < (i + 1); j++) {
            double sum = 0;
            for (int k = 0; k < j; k++)
                sum += L[i][k] * L[j][k];
            if (i == j)
                L[i][j] = sqrt(A[i][i] - sum);
            else
                L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
        }
    }
    return L;
}

vec Cholesky::solveTriangleSLE(const Matrix &A, const vec &B, MATRIX_TYPE type) {
    int dim = B.size();
    vec X(dim);
    if (type == MATRIX_TYPE::LOWER_TRIANGLE) {
        for (int i = 0; i < dim; ++i) {
            double sum = 0;
            for (int k = 0; k < i; ++k)
                sum += A[i][k] * X[k];
            X[i] = (B[i] - sum) / A[i][i];
        }
    } else if (type == MATRIX_TYPE::UPPER_TRIANGLE) {
        for (int i = dim - 1; i > -1; --i) {
            double sum = 0;
            for (int k = i + 1; k < dim; ++k)
                sum += A[i][k] * X[k];
            X[i] = (B[i] - sum) / A[i][i];
        }
    }
    return X;
}

vec Cholesky::CholeskySLE(const Matrix &A, const vec &B) {
    Matrix L = CholeskyDecomp(A);
    vec y = solveTriangleSLE(L, B, MATRIX_TYPE::LOWER_TRIANGLE);
    return solveTriangleSLE(transpose(L), y, MATRIX_TYPE::UPPER_TRIANGLE);
}

