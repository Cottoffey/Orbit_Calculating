#pragma once

#include <iostream>
#include <vector>


using vec = std::vector<double>;
using Matrix = std::vector<vec>;
enum class MATRIX_TYPE {
    LOWER_TRIANGLE = 0,
    UPPER_TRIANGLE,
};


void printMatrix(const Matrix &A);

Matrix transpose(Matrix &A);

vec multMatrix(const Matrix &A, const vec &B);

Matrix multMatrix(const Matrix &A, const Matrix &B);

