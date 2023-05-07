#include "Matrix.h"

void printMatrix(const Matrix &A) {
    for (const auto &y: A) {
        for (double x: y)
            std::cout << x << ' ';
        std::cout << '\n';
    }
}

Matrix transpose(Matrix &A) {
    int n = A.size();
    int m = A[0].size();
    Matrix C(m);
    for (int i = 0; i < m; ++i)
        C[i].resize(n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j)
            C[j][i] = A[i][j];
    }
    return C;

}


vec multMatrix(const Matrix &A, const vec &B) {
    int n = A.size();
    int m = B.size();
    vec C(n);

    for (size_t y = 0; y < n; ++y) {
        for (size_t x = 0; x < m; ++x)
            C[y] += A[y][x] * B[x];
    }
    return C;
}

Matrix multMatrix(const Matrix &A, const Matrix &B) {
    int n = A.size();
    int m1 = B.size();
    int m2 = B[0].size();
    Matrix C(n);

    for (int y = 0; y < n; ++y) {
        C[y].resize(m2);
        for (int x = 0; x < m2; ++x) {
            for (int k = 0; k < m1; ++k)
                C[y][x] += A[y][k] * B[k][x];
        }
    }
    return C;
}
