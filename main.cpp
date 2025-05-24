#include <omp.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "matrix.h"

int main() {
    omp_set_num_threads(4);  // Фиксируем число потоков
    omp_set_dynamic(0);      // Запрещаем автоматическое изменение числа потоков


    std::string filename = "bcsstk16.mtx";
    /*
    не симметричная      tols90.mtx
    симметричные         bcspwr01.mtx
    нет диагонали        bcsstk08.mtx  bcsstk11.mtx
    симметричные полож.  bcsstm08.mtx  bcsstk13.mtx  bcsstk16.mtx  bcsstk01.mtx  bcsstk04.mtx
    */
  /*
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

    //std::cout << "-----------------------------------------------\n";

    CSRMatrix A = convert_to_csr(mat);
    //print_csr_matrix(A);
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
    std::cout << "Some x, that was found: " << std::fixed << std::setprecision(10) << x[0] << ", " << x[1] << ", " << x[A.n_rows-1] << std::endl;

    // std::cout << "[";
    // for (size_t i = 0; i < x.size(); ++i) {
    //     std::cout << std::fixed << std::setprecision(10) << x[i];
    //     if (i != x.size() - 1) std::cout << ", ";
    // }
    // std::cout << "]" << std::endl;

    double E = calculate_mean_error(x, n, A, b);
    //std::cout << "\nAverage error: " << E << std::endl;

 */
 // /*
    double start, end;
    double n = 6.000;
    
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

        // Замер generate_b
    start = omp_get_wtime();
    std::vector<double> gb = generate_b(csr, n);
    end = omp_get_wtime();
    printf("generate_b: %.3f ms\n", (end-start)*1000);

        // Замер solve_ldlt
    start = omp_get_wtime();
    std::vector<double> x = solve_ldlt(csr, gb);
    end = omp_get_wtime();
    printf("solve_ldlt: %.3f ms\n", (end-start)*1000);

    // Замер calculate_mean_error
    start = omp_get_wtime();
    double cme = calculate_mean_error(x, n, csr, gb);
    end = omp_get_wtime();
    printf("calculate_mean_error: %.3f ms\n", (end-start)*1000);
 // */

    return 0;
}


// g++ -mconsole -o my_program.exe *.cpp
// ./my_program.exe
// g++ -fopenmp -o my_program.exe main.cpp matrix.cpp