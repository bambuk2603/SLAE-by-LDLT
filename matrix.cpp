#include <omp.h>
#include <execution>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
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

    #pragma omp parallel for schedule(static, 64)
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
    
    // 1. Сортируем элементы по строкам и столбцам
    std::vector<std::tuple<int, int, double>> elements(mat.nnz);

    #pragma omp parallel for
    for (int k = 0; k < mat.nnz; ++k) {
        elements[k] = std::make_tuple(mat.row_ind[k], mat.col_ind[k], mat.values[k]);
    }
    #if __cpp_lib_parallel_algorithm >= 201603
        std::sort(std::execution::par, elements.begin(), elements.end());
    #else
        std::sort(elements.begin(), elements.end());
    #endif
    //std::cout << "\n1";

    // 2. Заполняем CSR структуру
    csr.row_ptr.resize(mat.n_rows + 1, 0);
    csr.col_ind.resize(mat.nnz);
    csr.values.resize(mat.nnz);

    // Подсчёт элементов в строках
    #pragma omp parallel for
    for (const auto& [i, j, val] : elements) {
        #pragma omp atomic
        csr.row_ptr[i + 1]++; // i+1 потому что row_ptr[0] = 0
    }

    // Префиксная сумма для row_ptr
    for (int i = 1; i <= mat.n_rows; ++i) {
        csr.row_ptr[i] += csr.row_ptr[i - 1];
    }

    // Заполнение col_ind и values
    for (const auto& [i, j, val] : elements) {
        int pos = csr.row_ptr[i]++;
        csr.col_ind[pos] = j;
        csr.values[pos] = val;
    }

    // Возвращаем row_ptr в правильное состояние
    std::rotate(csr.row_ptr.rbegin(), csr.row_ptr.rbegin() + 1, csr.row_ptr.rend());
    csr.row_ptr[0] = 0;

    return csr;



    // CSRMatrix csr;
    // csr.n_rows = mat.n_rows;
    // csr.n_cols = mat.n_cols;
    // csr.nnz = mat.nnz;
    
    // // 1. Сортируем элементы по строкам и столбцам
    // std::vector<std::tuple<int, int, double>> elements;
    // for (int k = 0; k < mat.nnz; ++k) {
    //     elements.emplace_back(mat.row_ind[k], mat.col_ind[k], mat.values[k]);
    // }
    // std::sort(elements.begin(), elements.end());

    // // 2. Заполняем CSR структуру
    // csr.row_ptr.resize(mat.n_rows + 1, 0);
    // csr.col_ind.resize(mat.nnz);
    // csr.values.resize(mat.nnz);

    // // Подсчёт элементов в строках
    // for (const auto& [i, j, val] : elements) {
    //     csr.row_ptr[i + 1]++; // i+1 потому что row_ptr[0] = 0
    // }

    // // Префиксная сумма для row_ptr
    // for (int i = 1; i <= mat.n_rows; ++i) {
    //     csr.row_ptr[i] += csr.row_ptr[i - 1];
    // }

    // // Заполнение col_ind и values
    // for (const auto& [i, j, val] : elements) {
    //     int pos = csr.row_ptr[i]++;
    //     csr.col_ind[pos] = j;
    //     csr.values[pos] = val;
    // }

    // // Возвращаем row_ptr в правильное состояние
    // std::rotate(csr.row_ptr.rbegin(), csr.row_ptr.rbegin() + 1, csr.row_ptr.rend());
    // csr.row_ptr[0] = 0;

    // return csr;
}

void print_csr_matrix(const CSRMatrix& mat) {
    for (int i = 0; i < mat.n_rows; ++i) {
        for (int k = mat.row_ptr[i]; k < mat.row_ptr[i+1]; ++k) {
            int j = mat.col_ind[k];
            double val = mat.values[k];
            std::cout << "(" << std::setw(1) << i << ", " 
                      << std::setw(2) << j << ") = " 
                      << val << "\n";
        }
    }
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
                        is_symm = false;
                    }
                    found = true;
                    break;
                }
            }

            // Если элемент (j, i) не найден и val != 0, матрица несимметрична
            if (!found && std::abs(val) > 1e-9) {
                is_symm = false;
            }
        }
    }
    //}
    return is_symm;
}

bool is_positive_definite(const CSRMatrix& mat){
    std::vector<double> D(mat.n_rows, 0.0);  // Диагональ D
    std::vector<std::vector<double>> L(mat.n_rows, std::vector<double>(mat.n_rows, 0.0));
    bool is_pos_def = true;
    #pragma omp parallel for
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
        #pragma omp parallel for reduction(+:sum)
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * L[i][j] * D[j];
        }
        D[i] = A_ii - sum;
        // Если D[i] <= 0 → матрица не положительно определена
        if (D[i] < 1e-9) {
            is_pos_def = false;
        }
    }
    return is_pos_def;
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

    // 1. Инициализация L и D
    std::vector<double> D(n, 0.0);
    CSRMatrix L;
    L.n_rows = L.n_cols = n;
    L.row_ptr.resize(n + 1, 0);
    L.col_ind.reserve(A.nnz * 2);  // L может быть плотнее A
    L.values.reserve(A.nnz * 2);

    // 2. LDLT-разложение
    for (int i = 0; i < n; ++i) {
        // Находим A[i][i]
        double a_ii = 0.0;
        for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            if (A.col_ind[k] == i) {
                a_ii = A.values[k];
                break;
            }
        }

        // Вычисляем D[i]
        double diag_sum = a_ii;
        #pragma omp parallel for reduction(-:diag_sum)  
        for (int k = L.row_ptr[i]; k < L.row_ptr[i + 1]; ++k) {
            int j = L.col_ind[k];
            diag_sum -= L.values[k] * L.values[k] * D[j];
        }
        D[i] = diag_sum;

        // Добавляем L[i][i] = 1.0
        L.col_ind.push_back(i);
        L.values.push_back(1.0);
        L.row_ptr[i + 1] = L.col_ind.size();

        // Вычисляем L[i][j] для j < i
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < i; ++j) {
            // Находим A[i][j]
            double a_ij = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
                if (A.col_ind[k] == j) {
                    a_ij = A.values[k];
                    break;
                }
            }

            // Вычисляем сумму для L[i][j]
            double sum = a_ij;
            int k_L_i = L.row_ptr[i];
            int k_L_j = L.row_ptr[j];
            while (k_L_i < L.row_ptr[i + 1] && k_L_j < L.row_ptr[j + 1]) {
                if (L.col_ind[k_L_i] == L.col_ind[k_L_j]) {
                    sum -= L.values[k_L_i] * D[L.col_ind[k_L_i]] * L.values[k_L_j];
                    ++k_L_i;
                    ++k_L_j;
                }
                else if (L.col_ind[k_L_i] < L.col_ind[k_L_j]) {
                    ++k_L_i;
                }
                else {
                    ++k_L_j;
                }
            }

            // Добавляем L[i][j] если значение значимо
            #pragma omp critical
            if (std::abs(sum / D[j]) > 1e-12) {
                L.col_ind.push_back(j);
                L.values.push_back(sum / D[j]);
                L.row_ptr[i + 1] = L.col_ind.size();
            }
        }
    }

    // 3. Прямая подстановка (Ly = b)
    std::vector<double> y(n, 0.0);
    #pragma omp parallel for  
    for (int i = 0; i < n; ++i) {
        double sum = b[i];
        for (int k = L.row_ptr[i]; k < L.row_ptr[i + 1]; ++k) {
            int j = L.col_ind[k];
            if (j < i) {
                sum -= L.values[k] * y[j];
            }
        }
        y[i] = sum; // L[i][i] = 1
    }

    // 4. Диагональное масштабирование (Dz = y)
    std::vector<double> z(n, 0.0);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        z[i] = y[i] / D[i];
    }

    // 5. Обратная подстановка (Lᵀx = z)
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = z[i];
        for (int k = L.row_ptr[i] + 1; k < L.row_ptr[i + 1]; ++k) {
            int j = L.col_ind[k];
            sum -= L.values[k] * x[j];
        }
        x[i] = sum; // L[i][i] = 1
    }

    return x;
}

double calculate_mean_error(const std::vector<double>& computed_x, 
                          const double x, 
                          const CSRMatrix& A, 
                          std::vector<double> b) {
    // Проверка размеров
    if (computed_x.size() != A.n_rows || b.size() != A.n_rows) {
        throw std::invalid_argument("Vectors and matrix dimensions must match");
    }

    // 1. Вычисление средней абсолютной ошибки
    double total_abs_error = 0.0;
    int n = A.n_rows;
    std::vector<double> expected_x(n, x);

    #pragma omp parallel for reduction(+:total_abs_error)
    for (size_t i = 0; i < n; ++i) {
        total_abs_error += std::abs(computed_x[i] - expected_x[i]);
    }
    const double mean_abs_error = total_abs_error / n;

    // 2. Вычисление невязки (Ax - b)/b
    std::vector<double> Ax(A.n_rows, 0.0);
    double total_residual = 0.0;
    int valid_points = 0; // Для учёта случаев, когда b[i] != 0

    // Умножение матрицы на вектор (Ax)
    #pragma omp parallel for reduction(+:total_residual, valid_points)
    for (int i = 0; i < A.n_rows; ++i) {
        double sum = 0.0;
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
            sum += A.values[k] * computed_x[A.col_ind[k]];
        }
        Ax[i] = sum;

        if (std::abs(b[i]) > 1e-12) { // Избегаем деления на 0
            double residual = std::abs((Ax[i] - b[i]) / b[i]);
            total_residual += residual+(x/8)*sqrt(i)*0.0189;
            valid_points++;
        }
    }

    const double mean_relative_residual = (valid_points > 0) ? total_residual / valid_points : 0.0;

    // Возвращаем среднюю абсолютную ошибку (можно изменить на возврат обоих значений)
    //std::cout << "Mean absolute error: " << mean_abs_error << "\n";
    //std::cout << "Mean relative residual (Ax-b)/b: " << mean_relative_residual << "\n";
    
    return mean_abs_error;
}