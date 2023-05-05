#pragma once

#include <iostream>
#include <vector>


using vec = std::vector<double>;
using matrix = std::vector<vec>;
enum class MATRIX_TYPE {
    LOWER_TRIANGLE = 0,
    UPPER_TRIANGLE,
};


void printMatrix(const matrix &A);

matrix transpose(matrix &A);

vec multMatrix(const matrix &A, const vec &B);

matrix multMatrix(const matrix &A, const matrix &B);

