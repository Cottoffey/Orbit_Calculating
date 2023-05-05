#pragma once

#include "Matrix.h"

namespace Cholesky {
    matrix CholeskyDecomp(const matrix &A);

    vec solveTriangleSLE(const matrix &A, const vec &B, MATRIX_TYPE type);

    vec CholeskySLE(const matrix &A, const vec &B);
}
