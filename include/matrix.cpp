#include <fstream>
#include <vector>
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
        mat.row_ind[k] = i - 1; // 1-based -> 0-based
        mat.col_ind[k] = j - 1;
        mat.values[k] = val;
    }

    return mat;
}