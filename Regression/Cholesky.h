#pragma once

#include "Matrix.h"

namespace Cholesky {
    Matrix CholeskyDecomp(const Matrix &A);

    vec solveTriangleSLE(const Matrix &A, const vec &B, MATRIX_TYPE type);

    vec CholeskySLE(const Matrix &A, const vec &B);
}
