#include "algebra.h"
#include <random>
#include <stdexcept>
#include <optional>
#include <type_traits>
#include <format>
#include <iostream>

// Creation of namespace
namespace algebra
{
    // "============================================="
    // "              create_matrix                  "
    // "============================================="
    // function template for matrix creation
    template <typename T>
    MATRIX<T> create_matrix(std::size_t rows,
                            std::size_t columns,
                            std::optional<MatrixType> type,
                            std::optional<T> lowerBound,
                            std::optional<T> upperBound)
    {
        MATRIX<T> mat(rows, std::vector<T>(columns));
        if (rows == 0 || columns == 0){
            throw std::invalid_argument("number of rows and columns can not be zero");
        }
        switch (type.value_or(MatrixType::Zeros)) {
            case MatrixType::Zeros:
                for (int i = 0; i < (int)rows; ++i) {
                    for (int j = 0; j < (int)columns; ++j) {
                        mat[i][j] = 0;
                    }
                }
                break;
            case MatrixType::Ones:
                for (int i = 0; i < (int)rows; ++i) {
                    for (int j = 0; j < (int)columns; ++j) {
                        mat[i][j] = 1;
                    }
                }
                break;
            case MatrixType::Identity:
                if (rows != columns){
                    throw std::invalid_argument("Identity matrix should have equal number of rows and columns");
                }
                for (int i = 0; i < (int)rows; ++i) {
                    mat[i][i] = 1;
                }
                break;
            case MatrixType::Random:
                if (!lowerBound || !upperBound) {
                    throw std::invalid_argument("Random matrix requires lower and upper bounds");
                }
                if (*lowerBound >= *upperBound) {
                    throw std::invalid_argument("Upper bound must be greater than lower bound");
                }

                std::random_device rd;
                std::mt19937 gen(rd());
                if constexpr (std::is_integral_v<T>) {
                    std::uniform_int_distribution<T> dist(*lowerBound, *upperBound);
                    for (int i = 0; i < (int)rows; ++i) {
                        for (int j = 0; j < (int)columns; ++j) {
                            mat[i][j] = dist(gen);
                        }
                    }
                } else {
                    std::uniform_real_distribution<T> dist(*lowerBound, *upperBound);
                    for (int i = 0; i < (int)rows; ++i) {
                        for (int j = 0; j < (int)columns; ++j) {
                            mat[i][j] = dist(gen);
                        }
                    }
                }
                break;
        }
        return mat;
    }
    // Explicit instantiation for required types
    template MATRIX<int> create_matrix<int>(std::size_t rows, std::size_t columns, std::optional<MatrixType> type,
                                            std::optional<int> lowerBound, std::optional<int> upperBound);

    template MATRIX<double> create_matrix<double>(std::size_t rows, std::size_t columns, std::optional<MatrixType> type,
                                                  std::optional<double> lowerBound, std::optional<double> upperBound);

    template MATRIX<float> create_matrix<float>(std::size_t rows, std::size_t columns, std::optional<MatrixType> type,
                                                std::optional<float> lowerBound, std::optional<float> upperBound);



    template<typename T>
    void display(const MATRIX<T>& matrix){
        for (int i = 0; i < (int)(matrix.size()); ++i) {
            for (int j = 0; j < (int)(matrix[0].size()); ++j) {
                std::cout << std::format("{:^7}", matrix[i][j]) << "|";
            }
            std::cout << "\n";
        }
    }

    template<typename T>
    MATRIX<T> sum_sub(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB, std::optional<std::string> operation){
        if (matrixA.empty() && matrixB.empty()) {
            return MATRIX<T>{}; // Return empty matrix
        }

        // 2. Check if either matrix is empty (but not both)
        if (matrixA.empty() || matrixB.empty()) {
            throw std::invalid_argument("One matrix is empty while the other is not");
        }
        if ((matrixA.size() != matrixB.size()) || (matrixA[0].size() != matrixB[0].size())){
            throw std::invalid_argument("rows and columns for two matrix must be equal or can not be zero");
        }
        MATRIX<T> matrixC(matrixA.size(), std::vector<T>(matrixA[0].size()));
        if(operation == "sub"){
            for (int i = 0; i < (int)(matrixA.size()); ++i) {
                for (int j = 0; j < (int)(matrixA[0].size()); ++j) {
                   matrixC[i][j] = matrixA[i][j]- matrixB[i][j];
                }
            }
        }else{
            for (int i = 0; i < (int)(matrixA.size()); ++i) {
                for (int j = 0; j < (int)(matrixA[0].size()); ++j) {
                    matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
                }
            }
        }
        return matrixC;
    }
    template MATRIX<int> sum_sub<int>(const MATRIX<int>& matrixA, const MATRIX<int>& matrixB, std::optional<std::string> operation);
    template MATRIX<float> sum_sub<float>(const MATRIX<float>& matrixA, const MATRIX<float>& matrixB, std::optional<std::string> operation);
    template MATRIX<double> sum_sub<double>(const MATRIX<double>& matrixA, const MATRIX<double>& matrixB, std::optional<std::string> operation);


    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrix, T scalar){
        if (matrix.empty()) {
            return matrix;
//            throw std::invalid_argument("number of rows and columns can not be zero");
        }
        MATRIX<T> matrixC(matrix[0].size(), std::vector<T>(matrix.size()));
        for (int i = 0; i < (int)matrix.size(); ++i) {
            for (int j = 0; j < (int)matrix[0].size(); ++j) {
                matrixC[i][j] = matrix[i][j] * scalar;
            }
        }
        return matrixC;
    }
    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB){
        if (matrixA.empty() || matrixB.empty()) {
            throw std::invalid_argument("number of rows and columns can not be zero");
        }
        if (matrixA[0].size() != matrixB.size()){
            throw std::invalid_argument("number of rows and columns must match for multiply");
        }else{

            std::size_t rows = matrixA.size();
            std::size_t cols = matrixB[0].size();
            std::size_t inner = matrixA[0].size();

            MATRIX<T> result(rows, std::vector<T>(cols, T()));
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    for (std::size_t k = 0; k < inner; ++k) {
                        result[i][j] += matrixA[i][k] * matrixB[k][j];
                    }
                }
            }
            return result;
        }
    }
    template MATRIX<int> multiply(const MATRIX<int>&, int);
    template MATRIX<float> multiply(const MATRIX<float>&, float);
    template MATRIX<double> multiply(const MATRIX<double>&, double);

    template MATRIX<int> multiply(const MATRIX<int>&, const MATRIX<int>&);
    template MATRIX<float> multiply(const MATRIX<float>&, const MATRIX<float>&);
    template MATRIX<double> multiply(const MATRIX<double>&, const MATRIX<double>&);


    template<typename T>
    MATRIX<T> hadamard_product(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB){
        if (matrixA.empty() || matrixB.empty()){
            return matrixA;
        }
        if ((matrixA.size() != matrixB.size()) || (matrixA[0].size() != matrixB[0].size())){
            throw std::invalid_argument("number of rows and columns must match");
        }
        MATRIX<T> matrixC(matrixA.size(), std::vector<T>(matrixA[0].size()));
        for (int i = 0; i < (int)matrixA.size(); ++i) {
            for (int j = 0; j < (int)matrixA[0].size(); ++j) {
                matrixC[i][j] = matrixA[i][j] * matrixB[i][j];
            }
        }
        return matrixC;
    }
    template MATRIX<int> hadamard_product(const MATRIX<int>& matrixA, const MATRIX<int>& matrixB);

    template MATRIX<float> hadamard_product(const MATRIX<float>& matrixA, const MATRIX<float>& matrixB);
    template MATRIX<double> hadamard_product(const MATRIX<double>& matrixA, const MATRIX<double>& matrixB);


    template<typename T>
    MATRIX<T> transpose(const MATRIX<T>& matrix){
        if (matrix.empty()){
            return matrix;
        }
        MATRIX<T> matrixC(matrix[0].size(), std::vector<T>(matrix.size()));
        for (int i = 0; i < (int)matrix.size(); ++i) {
            for (int j = 0; j < (int)matrix[0].size(); ++j) {
                matrixC[j][i] = matrix[i][j];
            }
        }
        return matrixC;
    }
    template MATRIX<int> transpose(const MATRIX<int>& matrix);
    template MATRIX<float> transpose(const MATRIX<float>& matrix);
    template MATRIX<double> transpose(const MATRIX<double>& matrix);

    template<typename T>
    T trace(const MATRIX<T>& matrix){
        if (matrix.empty()){
            throw std::invalid_argument("must be square or can not be zero");
        }
        if (matrix.size() != matrix[0].size()){
            throw std::invalid_argument("must be square or can not be zero");
        }
        T result = 0;
        for (int i = 0; i < (int)matrix.size(); ++i) {
            result += matrix[i][i];
        }
        return result;
    }
    template int trace(const MATRIX<int>& matrix);
    template float trace(const MATRIX<float>& matrix);
    template double trace(const MATRIX<double>& matrix);


    template<typename T>
    T determinant(const MATRIX<T>& matrix) {
        int n = matrix.size();
        if (n == 0 || matrix[0].size() != static_cast<std::size_t>(n))
        {
            throw std::invalid_argument("Matrix must be square");
        }
        MATRIX<T> A = matrix;
        T det = 1;
        int swap_count = 0;

        for (int i = 0; i < n; ++i) {
            // Find pivot
            int pivot = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(A[j][i]) > std::abs(A[pivot][i])) {
                    pivot = j;
                }
            }
            // Swap rows if needed
            if (pivot != i) {
                std::swap(A[i], A[pivot]);
                swap_count++;
            }
            // If diagonal element is zero, det = 0
            if (A[i][i] == 0) {
                return 0;
            }
            det *= A[i][i];
            // Eliminate below
            for (int j = i + 1; j < n; ++j) {
                T factor = A[j][i] / A[i][i];
                for (int k = i; k < n; ++k) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }
        // Apply row swap effect
        if (swap_count % 2 != 0) {
            det *= -1;
        }
        return det;
    }

    template int determinant(const MATRIX<int>& matrix);
    template float determinant(const MATRIX<float>& matrix);
    template double determinant(const MATRIX<double>& matrix);

    template<typename T>
    MATRIX<T> inverse(const MATRIX<T>& input) {
        std::size_t n = input.size();
        if (n == 0 || input[0].size() != n) {
            throw std::invalid_argument("Matrix must be square and non-empty");
        }
        MATRIX<T> A = input;
        MATRIX<T> I(n, std::vector<T>(n, 0));
        for (std::size_t i = 0; i < n; ++i) {
            I[i][i] = 1;
        }
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t pivot = i;
            for (std::size_t j = i + 1; j < n; ++j) {
                if (std::abs(A[j][i]) > std::abs(A[pivot][i])) {
                    pivot = j;
                }
            }
            if (A[pivot][i] == 0) {
                throw std::runtime_error("Matrix is singular and cannot be inverted");
            }
            std::swap(A[i], A[pivot]);
            std::swap(I[i], I[pivot]);

            T div = A[i][i];
            for (std::size_t j = 0; j < n; ++j) {
                A[i][j] /= div;
                I[i][j] /= div;
            }

            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) {
                    T factor = A[j][i];
                    for (std::size_t k = 0; k < n; ++k) {
                        A[j][k] -= factor * A[i][k];
                        I[j][k] -= factor * I[i][k];
                    }
                }
            }
        }

        return I;
    }

    template MATRIX<int>inverse(const MATRIX<int>& input);
    template MATRIX<float>inverse(const MATRIX<float>& input);
    template MATRIX<double>inverse(const MATRIX<double>& input);
}
