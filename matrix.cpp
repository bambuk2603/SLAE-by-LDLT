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

    for (int k = 0; k < mat.nnz; ++k) {
        int i, j;
        double val;
        file >> i >> j >> val;
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
        
    if (line.find("symmetric") != std::string::npos) {
        return true;
    }
    else {
    // Быстрая проверка: матрица должна быть квадратной
        if (mat.n_rows != mat.n_cols) {
            return false;
        }
        else{
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
                                return false;
                            }
                            found = true;
                            break;
                        }
                    }

                    if (!found && val != 0.0) {  // Нулевые элементы могут быть пропущены
                        return false;
                    }
                }
            }
        }
    }
    return true;
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
        if (D[i] <= 1e-9) {
            return false;
        }
        else {
            return true;
        }
    }
}
