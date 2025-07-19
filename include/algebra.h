#ifndef AUT_AP_2024_Spring_HW1
#define AUT_AP_2024_Spring_HW1

#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <iostream>
#include <vector>
#include <optional>
#include <random>
#include "algebra.h"
#include <format>

namespace algebra
{
    template<typename T>
    using MATRIX = std::vector<std::vector<T>>;

// Matrix initialization types
    enum class MatrixType { Zeros, Ones, Identity, Random };
    // Function template for matrix initialization
    template<typename T>
    MATRIX<T> create_matrix(std::size_t rows,
                            std::size_t columns,
                            std::optional<MatrixType> type = MatrixType::Zeros,
                            std::optional<T> lowerBound = std::nullopt,
                            std::optional<T> upperBound = std::nullopt);

    template<typename T>
    void display(const MATRIX<T>& matrix);

    template<typename T>
    MATRIX<T> sum_sub(const MATRIX<T>& matrixA,
                      const MATRIX<T>& matrixB,
                      std::optional<std::string> operation = "sum");

    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrix, T scalar);
    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB);

    template<typename T>
    MATRIX<T> hadamard_product(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB);

    template<typename T>
    MATRIX<T> transpose(const MATRIX<T>& matrix);

    template<typename T>
    T trace(const MATRIX<T>& matrix);

    template<typename T>
    T determinant(const MATRIX<T>& matrix);

    template<typename T>
    MATRIX<T> inverse(const MATRIX<T>& input);
}
#endif
#endif