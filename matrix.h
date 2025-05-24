#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>
#include <vector>
#include <string>

struct Matrix {
    int n_rows, n_cols, nnz;
    std::vector<int> row_ind, col_ind;
    std::vector<double> values;
};

struct CSRMatrix {
    int n_rows, n_cols, nnz;
    std::vector<double> values;
    std::vector<int> col_ind, row_ptr;
};

Matrix read_mtx(const std::string& filename);
CSRMatrix convert_to_csr(const Matrix& mat);
void print_csr_matrix(const CSRMatrix& mat);
bool is_symmetric(const std::string& filename, const CSRMatrix& mat);
bool is_positive_definite(const CSRMatrix& mat);
std::vector<double> generate_b(const CSRMatrix& A, double n);
std::vector<double> solve_ldlt(const CSRMatrix& A, const std::vector<double>& b);
double calculate_mean_error(const std::vector<double>& computed_x, 
                          const double x, const CSRMatrix& A, std::vector<double> b);

#endif