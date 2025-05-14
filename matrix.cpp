#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

Matrix read_mtx(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("File not found");

    Matrix mat;
    std::string line;

    // Пропуск заголовка
    do { std::getline(file, line); } 
    while (line[0] == '%');

    // Чтение размеров
    sscanf(line.c_str(), "%d %d %d", &mat.n_rows, &mat.n_cols, &mat.nnz);

    // Чтение данных
    mat.row_ind.resize(mat.nnz);
    mat.col_ind.resize(mat.nnz);
    mat.values.resize(mat.nnz);
    
    std::vector<std::string> file_lines(mat.nnz);
    for (int k = 0; k < mat.nnz; ++k) {
        std::getline(file, file_lines[k]);
    }

    #pragma omp parallel for schedule(static, 1024)
    for (int k = 0; k < mat.nnz; ++k) {
        int i, j;
        double val;
        sscanf(file_lines[k].c_str(), "%d %d %lf", &i, &j, &val);
        mat.row_ind[k] = i-1; // 1-based -> 0-based
        mat.col_ind[k] = j-1;
        mat.values[k] = val;
    }

    return mat;
}

CSRMatrix convert_to_csr(const Matrix& mat) {
    CSRMatrix csr;
    csr.n_rows = mat.n_rows;
    csr.n_cols = mat.n_cols;
    csr.nnz = mat.nnz;
    csr.values = mat.values;
    csr.col_ind = mat.col_ind;
    csr.row_ptr.resize(mat.n_rows + 1, 0);

    // Подсчёт элементов в строках
    for (int k = 0; k < mat.nnz; ++k) 
        csr.row_ptr[mat.row_ind[k]]++;

    // Префиксная сумма для row_ptr
    std::partial_sum(csr.row_ptr.begin(), csr.row_ptr.end(), csr.row_ptr.begin());
    std::rotate(csr.row_ptr.rbegin(), csr.row_ptr.rbegin() + 1, csr.row_ptr.rend());
    csr.row_ptr[0] = 0;

    return csr;
}

bool is_symmetric(const std::string& filename, const CSRMatrix& mat) {
    // Сначала проверяем заголовок файла на наличие symmetric
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("File not found");

    std::string line;
    do { std::getline(file, line); }
        while (line[0] == '\n');

            // Быстрая проверка: матрица должна быть квадратной
        if (mat.n_rows != mat.n_cols) {
            return false;
        }

        if (line.find("symmetric") != std::string::npos) {
            return true;
        }

        bool is_symm = true;

        //#pragma omp parallel shared(is_symm)
        //{
        //#pragma omp for nowait
        // Проверяем симметричность данных в CSR формате
        for (int i = 0; i < mat.n_rows; ++i) {
            for (int k = mat.row_ptr[i]; k < mat.row_ptr[i+1]; ++k) {
                int j = mat.col_ind[k];
                double val = mat.values[k];
                // Ищем элемент A[j][i]
                bool found = false;
                for (int m = mat.row_ptr[j]; m < mat.row_ptr[j+1]; ++m) {
                    if (mat.col_ind[m] == i) {
                        if (std::abs(mat.values[m] - val) > 1e-9) {
                            //#pragma omp atomic write
                            is_symm = false;
                        }
                        found = true;
                        break;
                    }
                    if (!found && val != 0.0) {
                    //#pragma omp atomic write
                    is_symm = false;
                    found = false;
                    }
                }
            }
        }
    //}
    return is_symm;
}

bool is_positive_definite(const CSRMatrix& mat){
    std::vector<double> D(mat.n_rows, 0.0);  // Диагональ D
    std::vector<std::vector<double>> L(mat.n_rows, std::vector<double>(mat.n_rows, 0.0));
    for (int i = 0; i < mat.n_rows; ++i) {
        // Вычисляем L[i][k] и D[i]
        double sum = 0.0;
        // Ищем A[i][i] в CSR
        double A_ii = 0.0;
        for (int k = mat.row_ptr[i]; k < mat.row_ptr[i + 1]; ++k) {
            if (mat.col_ind[k] == i) {
                A_ii = mat.values[k];
                break;
            }
        }
        // Вычисляем сумму для D[i]
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * L[i][j] * D[j];
        }
        D[i] = A_ii - sum;
        // Если D[i] <= 0 → матрица не положительно определена
        if (D[i] > 1e-9) {
            return true;
        }
            return false;
    }
}

std::vector<double> generate_b(const CSRMatrix& A, double n) {
    std::vector<double> b(A.n_rows, 0.0);
    for (int i = 0; i < A.n_rows; ++i) {
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
            b[i] += A.values[k] * n;  // Все x[j] = n
        }
    }
    return b;
}

std::vector<double> solve_ldlt(const CSRMatrix& A, const std::vector<double>& b) {
    const int n = A.n_rows;

    for (int i = 0; i < A.n_rows; ++i) {
    bool has_diagonal = false;
    for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
        if (A.col_ind[k] == i) has_diagonal = true;
    }
    if (!has_diagonal) {
        std::cerr << "Error: Missing diagonal element for row " << i << std::endl;
        return {}; // Возвращаем пустой вектор
    }
}
    
    // 1. Преобразуем CSR в плотную нижнюю треугольную матрицу
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<double> D(n, 0.0);
    
    // Заполняем нижний треугольник из CSR
    for (int i = 0; i < n; ++i) {
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
            int j = A.col_ind[k];
            if (j <= i) {  // Только нижний треугольник
                L[i][j] = A.values[k];
            }
        }
    }

    // 2. LDLᵀ-разложение (на месте)
    for (int j = 0; j < n; ++j) {
        // Вычисляем D[j]
        double sum = 0.0;
        for (int k = 0; k < j; ++k) {
            sum += L[j][k] * L[j][k] * D[k];
        }
        D[j] = L[j][j] - sum;
        L[j][j] = 1.0;  // Диагональ L

        // Вычисляем L[i][j] для i > j
        for (int i = j+1; i < n; ++i) {
            sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum += L[i][k] * L[j][k] * D[k];
            }
            L[i][j] = (L[i][j] - sum) / D[j];
        }
    }

    // 3. Решение Ly = b (прямая подстановка)
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int k = 0; k < i; ++k) {
            y[i] -= L[i][k] * y[k];
        }
    }

    // 4. Решение Dz = y
    std::vector<double> z(n);
    for (int i = 0; i < n; ++i) {
        z[i] = y[i] / D[i];
    }

    // 5. Решение Lᵀx = z (обратная подстановка)
    std::vector<double> x(n);
    for (int i = n-1; i >= 0; --i) {
        x[i] = z[i];
        for (int k = i+1; k < n; ++k) {
            x[i] -= L[k][i] * x[k];  // Используем Lᵀ
        }
    }

    return x;
}
