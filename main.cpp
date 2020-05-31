//
// Created by Ahtyamov Danil on 28.04.2020.
// https://github.com/ahtyamovdanil
//

#include <mpi.h>
#include <iostream>
#include <vector>
#include <array>
#include <cassert>
#include <cmath>
#include "solver.hpp"

#define N 6
#define M 4

//  mpiCC -o out main.cpp && mpiexec -np 8 /home/mirror/university/diploma/SimplexSolver/out

int main(int argc, char *argv[]) {
    int numprocs, myid, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    simplex::init_solver(&argc, &argv, &numprocs, &myid, &namelen, processor_name);

    //constexpr int N = 6;
    //constexpr int M = 4;

    typedef  std::array<std::array<float, N>, M> Matrix;
    //typedef  std::array<std::array<float, M>, N> MatrixT;
    //typedef  std::vector<std::vector<float>> Matrix;
    typedef  std::vector<float> MatrixT;

    // number of rows in one batch
    int batch_rows = std::ceil(float(N) / numprocs);
    // batch size
    int batch_size = batch_rows * M;
    // last batch may be smaller than the other ones
    int rows_in_last = (N % batch_rows == 0) ? batch_rows : N % batch_rows;

    // constraint matrix
    /*
    Matrix A = {{{6,4,1,0,0,0},
                 {1,2,0,1,0,0},
                 {-1,1,0,0,1,0},
                 {0,1,0,0,0,1}}};
*/
    // transposed constraint matrix
    /*
    MatrixT At = {{{6, 1, -1, 0},
                   {4, 2, 1, 1},
                   {1, 0, 0, 0},
                   {0, 1, 0, 0},
                   {0, 0, 1, 0},
                   {0, 0, 0, 1}}};
    */
    MatrixT At = {6, 1, -1, 0,
                  4, 2, 1, 1,
                  1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1};

    std::vector<float> C = {-5, -4, 0, 0, 0, 0}; // z row
    std::vector<float> B = {24, 6, 1, 2};    // F(x)
    float optimum = 0;
    int pivot_row_idx;
    std::vector<float> pivot_column;

    // Matrix with bunch of columns for each processor
    std::vector<std::array<float, M>> slave_columns(batch_rows);

    // last process may have smaller batch then other ones
    if (myid == 0) {
        slave_columns.resize(rows_in_last);
    }

    // Vector with local pivots from all processors
    std::vector<float> local_piv_values(numprocs);
    std::vector<int> local_piv_idxs(numprocs);

    // Vector with local z row
    std::vector<float> local_c(batch_rows);
    if (myid == 0) {
        local_c.resize(rows_in_last);
    }

    // indicator that optimum has been found
    bool is_done = false;

    //broadcast vector [B] to all processors
    MPI_Bcast(B.data(), M, MPI_FLOAT, 0, MPI_COMM_WORLD);

    assert(("Number of processors must be less or equal to the number of variables", N >= numprocs));
    //scatter rows of matrix [At] to each processor
    std::cout << "kek_0" << std::endl;
    MPI_Scatter(At.data(), batch_size, MPI_FLOAT, slave_columns.at(0).data(), batch_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //scatter elements of z-row [C] to each processor
    MPI_Scatter(C.data(), batch_rows, MPI_FLOAT, local_c.data(), batch_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    std::cout << "kek_1" << std::endl;
    std::vector<float>::iterator loc_piv_iter;
    int loc_piv_i;
    float loc_piv_val;

    std::vector<float>::iterator zvalue_iter;
    int pivot_column_proc_id;
    float global_zvalue;
    int global_zidx;

    while(!is_done){
        // find pivot column value and index in each processor
        loc_piv_iter = std::min_element(local_c.begin(), local_c.end());
        loc_piv_i = std::distance(local_c.begin(), loc_piv_iter);
        loc_piv_val = *loc_piv_iter;
        pivot_column = slave_columns.at(loc_piv_i);
        std::cout << "kek_2" << std::endl;
        // gather all local pivots values and indexes from all processors
        MPI_Gather(&loc_piv_val, 1, MPI_FLOAT, local_piv_values.data(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Gather(&loc_piv_i, 1, MPI_FLOAT, local_piv_idxs.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        std::cout << "kek_3" << std::endl;
        // get value of global pivot column and its processor id
        if (myid == 0) {
            zvalue_iter = std::min_element(local_piv_values.begin(), local_piv_values.end());
            pivot_column_proc_id = std::distance(local_piv_values.begin(), zvalue_iter);
            global_zvalue = *zvalue_iter;
            global_zidx = local_piv_idxs[pivot_column_proc_id];
        }

        // broadcast processor id to all processors
        MPI_Bcast(&pivot_column_proc_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "kek_4" << std::endl;
        // broadcast pivot column to all processors
        MPI_Bcast(pivot_column.data(), M, MPI_FLOAT, pivot_column_proc_id, MPI_COMM_WORLD);
        std::cout << "kek_5" << std::endl;
        // broadcast pivot value in the z-row
        MPI_Bcast(&global_zvalue, 1, MPI_INT, 0, MPI_COMM_WORLD);
        std::cout << "kek_6" << std::endl;

        if (global_zvalue >= 0) {
            is_done = true;
        }

        //find index of the pivot row
        if (myid == 0) {
            pivot_row_idx = simplex::find_piv_row(B, pivot_column);
            std::cout << "pivot_column" << std::endl;
            for(auto & item : pivot_column)
                std::cout << item << ' ';
            std::cout << std::endl;
        }
        std::cout << "kek_7" << std::endl;
        MPI_Bcast(&pivot_row_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
        std::cout << "kek_8" << std::endl;
        std::cout << pivot_row_idx << std::endl;
        // pivot value - element on the cross of the pivot row and pivot column
        float pivot_elem = pivot_column[pivot_row_idx];

        std::cout << "kek_10" << std::endl;
        // update pivot row in each processor
        for (int i=0; i<slave_columns.size(); ++i) {
            slave_columns[i][pivot_row_idx] /= pivot_elem;
        }

        // update [B] in pivot row
        B[pivot_row_idx] /= pivot_elem;

        // update global optimum value
        if (myid == 0) {
            //std::cout << "kek0 " << global_zvalue << ' ' << B[pivot_row_idx] << std::endl;
            optimum -= global_zvalue * B[pivot_row_idx];
        }

        // update z-row
        for (int j = 0; j < local_c.size(); ++j) {
            local_c[j] -= global_zvalue * slave_columns[j][pivot_row_idx];
        }

        //update values in the rest of the table
        for (int i=0; i<slave_columns.size(); ++i) {
            for (int j = 0; j < M; ++j) {
                if (j != pivot_row_idx) {
                    slave_columns[i][j] -= pivot_column[j] * slave_columns[i][pivot_row_idx];
                }
            }
        }

        // update results vector [B]
        for (int j = 0; j < M; ++j) {
            if (j != pivot_row_idx) {
                B[j] -= pivot_column[j] * B[pivot_row_idx];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // gather vector C from all processors
    MPI_Gather(local_c.data(), batch_rows, MPI_FLOAT, C.data(), batch_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // gather all columns
    MPI_Gather(slave_columns.data(), batch_size, MPI_FLOAT, At.data(), batch_size, MPI_INT, 0, MPI_COMM_WORLD);

    if(myid==0){
        //simplex::print_matrix<N, M>("results", At_test);
        /*
        std::cout << "local pivots values:" << std::endl;
        for(float &item: local_piv_values)
            std::cout << item << ' ';

        std::cout << "\nlocal pivots indexes:" << std::endl;
        for(int &item: local_piv_idxs)
            std::cout << item << ' ';

        std::cout << "\npivot column:" << std::endl;
        for(float &item: pivot_column)
            std::cout << item << ' ';

        //std::cout << std::endl << "global pivot:" << std::endl;
        std::cout << "\nglobal pivot value: " << global_zvalue << std::endl
                  << "global pivot processor id: " << pivot_column_proc_id << std::endl
                  << "global pivot column index: : " << global_zidx << std::endl
                  << "pivot row index: : " << pivot_row_idx << std::endl;
*/
        std::cout << "----------------------" <<std::endl;
        std::cout << "\nB:" << std::endl;

        for(int i=0; i<B.size(); ++i)
            std::cout << B[i] << ' ';
        std::cout << "\nC:" << std::endl;
        for(int i=0; i<C.size(); ++i)
            std::cout << C[i] << ' ';
        //simplex::print_matrix("\nresult matrix: ", At);


        std::cout << "optimum value: " << optimum << std::endl;
        std::cout << "ok" << std::endl;
    }

    // Освобождение подсистемы MPI
    MPI_Finalize();
}
