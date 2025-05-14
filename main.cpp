#include <omp.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "matrix.h"

int main() {
    std::string filename = "bcsstm08.mtx";
    /*
    не симметричная      tols90.mtx
    симметричные         bcspwr01.mtx
    симметричные полож.  bcsstk01.mtx  bcsstm08.mtx  bcsstk04.mtx
    нет диагонали        bcsstk08.mtx  bcsstk11.mtx
    подходят             bcsstm08.mtx
    */
// /*
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

    CSRMatrix A = convert_to_csr(mat);
/////////////////////////////////////////////////////////
    bool symmetric = is_symmetric(filename, A);
        if (symmetric) {std::cout << "\nMatrix is symmetric!";}
        else {std::cout << "\nMatrix is not symmetric :(";
              return 0; };

    bool positive_definite = is_positive_definite(A);
        if (positive_definite) {std::cout << "\nMatrix is positive definite!";}
        else {std::cout << "\nMatrix is not positive definite :(" << std::endl;
              return 0; };
    double n;
    std::cout << "\nEnter x: ";
    std::cin >> n;
    std::vector<double> b = generate_b(A, n);
    // std::cout << b[0] << ", " << b[1472] << std::endl;

    std::cout << "\nI've generated right part of SLAE. Lets solve!\n";

    std::vector<double> x = solve_ldlt(A, b);
    std::cout << "Some x, that was found: " << std::fixed << std::setprecision(10) << x[0] << ", " << x[A.n_rows-1] << std::endl;

// */
/*
    double start, end;
    
    // Замер read_mtx
    start = omp_get_wtime();
    Matrix mat = read_mtx(filename);
    end = omp_get_wtime();
    printf("read_mtx: %.3f ms\n", (end-start)*1000);
    
    // Замер convert_to_csr
    start = omp_get_wtime();
    CSRMatrix csr = convert_to_csr(mat);
    end = omp_get_wtime();
    printf("convert_to_csr: %.3f ms\n", (end-start)*1000);
    
    // Замер is_symmetric
    start = omp_get_wtime();
    bool symm = is_symmetric(filename, csr);
    end = omp_get_wtime();
    printf("is_symmetric: %.3f ms\n", (end-start)*1000);
    
    // Замер is_positive_definite
    start = omp_get_wtime();
    bool pd = is_positive_definite(csr);
    end = omp_get_wtime();
    printf("is_positive_definite: %.3f ms\n", (end-start)*1000);
*/

    return 0;
}


// g++ -mconsole -o my_program.exe *.cpp
// ./my_program.exe
// g++ -fopenmp -o my_program.exe main.cpp matrix.cpp