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


//#define N 6
//#define M 4

//  mpiCC -o out main.cpp && mpiexec -np 8 /home/mirror/university/diploma/SimplexSolver/out

int main(int argc, char *argv[]) {

    int numprocs, myid, namelen;

    char processor_name[MPI_MAX_PROCESSOR_NAME];

    solver::init_solver(&argc, &argv, &numprocs, &myid, &namelen, processor_name);

    int const N = 60000;
    int const M = 100;
    // transposed constraint matrix

    /*
    solver::matrix<float> At(N, M, {6, 1, -1, 0,
                                    4, 2, 1, 1,
                                    1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, 1});

    std::vector<float> C = {-5, -4, 0, 0, 0, 0}; // z row
    std::vector<float> B = {24, 6, 1, 2};    // F
*/
// todo
    solver::matrix<float> At(N, M);
    std::vector<float> B(M), C(N);

    if(myid == 0)
        std::tie(At, B, C) = solver::random_problem<float>(N, M);

// todo



    //int const& qwe = numprocs;
    //constexpr int N = 6;
    //constexpr int M = 4;

    //typedef  std::array<std::array<float, N>, M> Matrix;
    //typedef  std::array<std::array<float, M>, N> MatrixT;
    //typedef  std::vector<std::vector<float>> Matrix;
    //  std::vector<float> MatrixT;


    // number of rows in one batch
    int batch_rows = std::floor(float(N) / numprocs);

    // batch size
    int batch_size = batch_rows * M;
    int big_batch_num = N % numprocs;
    // last batch may be smaller than the other ones

    std::vector<int> batch_sizes(numprocs);
    std::vector<int> batch_offsets(numprocs);
    std::vector<int> c_sizes(numprocs);
    std::vector<int> c_offsets(numprocs);

    for(int i=0; i<numprocs; ++i){
        if(i>=numprocs-big_batch_num){
            batch_sizes[i] = batch_size+M;
            c_sizes[i] = batch_rows+1;
            batch_offsets[i] = batch_offsets[i-1] + batch_sizes[i-1];
            c_offsets[i] = c_offsets[i-1] + c_sizes[i-1];
        }
        else{
            batch_sizes[i] = batch_size;
            c_sizes[i] = batch_rows;
            if(i==0){
                batch_offsets[i] = 0;
                c_offsets[i] = 0;
            }
            else{
                batch_offsets[i] = batch_offsets[i-1] + batch_sizes[i-1];
                c_offsets[i] = c_offsets[i-1] + c_sizes[i-1];
            }
        }
    }
    // constraint matrix
    /*
    Matrix A = {{{6,4,1,0,0,0},
                 {1,2,0,1,0,0},
                 {-1,1,0,0,1,0},
                 {0,1,0,0,0,1}}};
*/


    float optimum = 0;
    int pivot_row_idx;
    std::vector<float> pivot_column;

    // Matrix with bunch of columns for each processor
    solver::matrix<float> slave_columns(batch_rows+1, M);

    // last process may have smaller batch then other ones
    slave_columns.resize(c_sizes[myid]);

    // Vector with local pivots from all processors
    std::vector<float> local_piv_values(numprocs);
    std::vector<int> local_piv_idxs(numprocs);

    // Vector with local z row
    std::vector<float> local_c(batch_rows+1);
    local_c.resize(c_sizes[myid]);

    // indicator that optimum has been found
    bool is_done = false;

    //broadcast vector [B] to all processors
    MPI_Bcast(B.data(), M, MPI_FLOAT, 0, MPI_COMM_WORLD);

    assert(("Number of processors must be less or equal to the number of variables", N >= numprocs));
    //scatter rows of matrix [At] to each processor
    MPI_Scatterv(At.data(), batch_sizes.data(), batch_offsets.data(), MPI_FLOAT, slave_columns.data(), batch_size+M, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(myid == 0){
        At.clear();
    }
    //scatter elements of z-row [C] to each processor
    MPI_Scatterv(C.data(), c_sizes.data(), c_offsets.data(), MPI_FLOAT, local_c.data(), batch_rows+1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    std::vector<float>::iterator loc_piv_iter;
    int loc_piv_i;
    float loc_piv_val;

    std::vector<float>::iterator zvalue_iter;
    int pivot_column_proc_id;
    float global_zvalue;
    int global_zidx;

    int i = 1;
    while(!is_done){
        // find pivot column value and index in each processor
        //for(auto& item: local_c)
        //    std::cout << item << ' ';

        loc_piv_iter = std::min_element(local_c.begin(), local_c.end());
        loc_piv_i = std::distance(local_c.begin(), loc_piv_iter);
        loc_piv_val = *loc_piv_iter;
        pivot_column = slave_columns.at(loc_piv_i);

        // gather all local pivots values and indexes from all processors
        MPI_Gather(&loc_piv_val, 1, MPI_FLOAT, local_piv_values.data(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Gather(&loc_piv_i, 1, MPI_FLOAT, local_piv_idxs.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        // get value of global pivot column and its processor id
        if (myid == 0) {
            // output for debuging
            //std::cout << "\n+-------------------------------+" << std::endl;
            std::cout << "|  current optimum: " << optimum << "\t\t|" << std::endl;
            /*
            std::cout << "|  loc mins: ";
            for(auto& item: local_piv_values)
                std::cout << item << ' ';
            */
            //std::cout << "\t\t|\n+-------------------------------+" << std::endl;

            zvalue_iter = std::min_element(local_piv_values.begin(), local_piv_values.end()-1);
            pivot_column_proc_id = std::distance(local_piv_values.begin(), zvalue_iter);
            global_zvalue = *zvalue_iter;
            global_zidx = local_piv_idxs[pivot_column_proc_id];
            ++i;

        }

        // broadcast processor id to all processors
        MPI_Bcast(&pivot_column_proc_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // broadcast pivot column to all processors
        MPI_Bcast(pivot_column.data(), M, MPI_FLOAT, pivot_column_proc_id, MPI_COMM_WORLD);

        // broadcast pivot value in the z-row
        MPI_Bcast(&global_zvalue, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (global_zvalue >= 0) {
            is_done = true;
        }

        //find index of the pivot row
        if (myid == 0) {
            pivot_row_idx = solver::find_piv_row(B, pivot_column);
        }

        MPI_Bcast(&pivot_row_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // pivot value = element on the cross of the pivot row and pivot column
        float pivot_elem = pivot_column[pivot_row_idx];
        // update pivot row in each processor
        for (int i=0; i<slave_columns.rows(); ++i) {
            slave_columns.at(i,pivot_row_idx) /= pivot_elem;
        }
        // update [B] in pivot row
        B[pivot_row_idx] /= pivot_elem;

        // update global optimum value
        if (myid == 0) {
            optimum -= global_zvalue * B[pivot_row_idx];
        }

        // update z-row
        for (int j = 0; j < local_c.size(); ++j) {
            auto abc = slave_columns.at(j,pivot_row_idx);
            local_c[j] -= global_zvalue * slave_columns.at(j,pivot_row_idx);
        }

        //update values in the rest of the table
        for (int i=0; i<slave_columns.rows(); ++i) {
            for (int j = 0; j < M; ++j) {
                if (j != pivot_row_idx) {
                    slave_columns.at(i,j) -= pivot_column[j] * slave_columns.at(i,pivot_row_idx);
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
    MPI_Gatherv(local_c.data(), batch_rows,  MPI_FLOAT, C.data(),c_sizes.data(), c_offsets.data(), MPI_FLOAT, 0, MPI_COMM_WORLD);

    // gather all columns
    //MPI_Gather(slave_columns.data(), batch_size, MPI_FLOAT, At.data(), batch_size, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // print finale result
    if(myid==0){
        std::cout << "----------------------" <<std::endl;
        std::cout << "\nB:" << std::endl;

        for(int i=0; i<B.size(); ++i)
            std::cout << B[i] << ' ';
        //std::cout << "\nC:" << std::endl;
        std::cout << "\noptimum value: " << optimum << std::endl;
        std::cout << "number of iterations: " << i << std::endl;
        std::cout << "done" << std::endl;
    }



    // Освобождение подсистемы MPI
    MPI_Finalize();
    return 0;
}
