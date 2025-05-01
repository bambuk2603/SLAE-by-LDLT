#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "matrix.h"

int main() {
    std::string filename = "bcspwr01.mtx";
    /*
    не симметричная   tols90.mtx
    симметричные      bcspwr01.mtx  bcsstk01.mtx
    */

    // Чтение матрицы
    Matrix mat = read_mtx(filename);

    // Вывод информации
    std::cout << "Rows: " << mat.n_rows 
              << ", Cols: " << mat.n_cols 
              << ", Non-zeros: " << mat.nnz << "\n";

    // Пример вывода первых 5 элементов
    std::cout << "First 5 values:\n";
    for (int i = 0; i < std::min(5, mat.nnz); ++i) {
        std::cout << "(" << mat.row_ind[i] << ", " 
                  << mat.col_ind[i] << ") = " 
                  << mat.values[i] << "\n";
    }

    CSRMatrix mat_csr = convert_to_csr(mat);
/////////////////////////////////////////////////////////
    bool symmetric;
    symmetric = is_symmetric(filename, mat_csr);
        if (symmetric) {std::cout << "\nMatrix is symmetric!";}
        else {std::cout << "\nMatrix is not symmetric :(";
              return 0; };

    bool positive_definite = is_positive_definite(mat_csr);
        if (positive_definite) {std::cout << "\nMatrix is positive definite!";}
        else {std::cout << "\nMatrix is not positive definite :(";
              return 0; };


    return 0;
}


// g++ -mconsole -o my_program.exe *.cpp
// ./my_program.exe