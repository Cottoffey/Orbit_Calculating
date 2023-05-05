#include "Matrix.h"

void printMatrix(const matrix &A) {
    for (const auto &y: A) {
        for (double x: y)
            std::cout << x << ' ';
        std::cout << '\n';
    }
}

matrix transpose(matrix &A) {
    matrix C;
    int n = A.size();
    int m = A[0].size();
    C.resize(m);
    for (int i = 0; i < m; ++i)
        C[i].resize(n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j)
            C[j][i] = A[i][j];
    }
    return C;

}


vec multMatrix(const matrix &A, const vec &B) {
    vec C;
    int n = A.size();
    int m = B.size();
    C.resize(n);

    for (size_t y = 0; y < n; ++y) {
        for (size_t x = 0; x < m; ++x)
            C[y] += A[y][x] * B[x];
    }
    return C;
}

matrix multMatrix(const matrix &A, const matrix &B) {
    matrix C;
    int n = A.size();
    int m1 = B.size();
    int m2 = B[0].size();
    C.resize(n);

    for (int y = 0; y < n; ++y) {
        C[y].resize(m2);
        for (int x = 0; x < m2; ++x) {
            for (int k = 0; k < m1; ++k)
                C[y][x] += A[y][k] * B[k][x];
        }
    }
    return C;
}
