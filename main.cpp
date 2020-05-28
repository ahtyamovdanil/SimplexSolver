#include <mpi.h>
#include <iostream>
#include <vector>
#include <array>
#include <cassert>
#include <cmath>
#include "solver.hpp"


//  mpiCC -o out main.cpp && mpiexec -np 8 /home/mirror/university/diploma/SimplexSolver/out

int main(int argc, char *argv[])
{
    int numprocs, myid, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    simplex::init_solver(&argc, &argv, &numprocs, &myid, &namelen, processor_name);

    //std::cout << "number of processors: " << numprocs << std::endl;
    constexpr int N = 5;
    constexpr int M = 4;

    using Matrix = std::array<std::array<float, N>, M>;
    using MatrixT = std::array<std::array<float, M>, N>;

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
    const MatrixT At = {{{6,1,-1,0},
                       {4,2,1,1},
                       {1,0,0,0},
                       {0,1,0,0},
                      // {0,0,1,0},
                       {0,0,0,1}}};

    MatrixT At_test;
    std::array<float, N> C = {-1,-1,1,-5,2}; // z row
    std::array<float, M> B = {24,6,1,2};      // F(x)

    // Matrix with bunch of columns for each processor
    std::vector<std::array<float, M>> slave_columns(batch_rows);
    // Vector with local pivots from all processors
    std::vector<float> local_pivots(numprocs);
    // Vector with local z row
    std::vector<float> local_c(batch_rows);


    //broadcast vector [C] to all processors
    MPI_Bcast(C.data(), N, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //broadcast vector [B] to all processors
    MPI_Bcast(B.data(), M, MPI_FLOAT, 0, MPI_COMM_WORLD);

    assert(("Number of processors must be less or equal to the number of variables", N >= numprocs));
    //scatter rows of matrix [At] to each processor
    MPI_Scatter(At.data(), batch_size, MPI_FLOAT, slave_columns.at(0).data(), batch_size, MPI_FLOAT, 0 , MPI_COMM_WORLD);

    // last process may have smaller batch then other ones
    if(myid == numprocs - 1){
        slave_columns.resize(rows_in_last);
    }

    //scatter elements of z-row [C] to each processor
    MPI_Scatter(C.data(), batch_rows, MPI_FLOAT, local_c.data(), batch_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if(myid == numprocs - 1){
        local_c.resize(rows_in_last);
    }

    // find pivot column in each processor
    float local_piv = *std::min_element(local_c.begin(), local_c.end());

    // gather all local pivots from all processors
    MPI_Gather(&local_piv, 1, MPI_FLOAT, local_pivots.data(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // get value of global pivot and its processor id
    auto pivot_iter = std::min_element(local_pivots.begin(), local_pivots.end());
    int pivot_id = std::distance(local_pivots.begin(), pivot_iter);
    float global_pivot = *pivot_iter;
    /*
    for(auto &col: slave_columns){
        for(auto &item: col) {
            item = myid;
            //std::cout << myid << ' ' << item << std::endl;
        }
    }
*/
    //MPI_Gather(slave_columns.at(0).data(), batch_size, MPI_INT, At_test.data(), batch_size, MPI_INT, 0, MPI_COMM_WORLD);

    if(myid==0){
        //simplex::print_matrix<N, M>("results", At_test);
        std::cout << "local pivots:" << std::endl;
        for(auto &item: local_pivots)
            std::cout << item << ' ';

        //std::cout << std::endl << "global pivot:" << std::endl;
        std::cout << "\nglobal pivot value: " << global_pivot << std::endl
                  << "global pivot processor id: " << pivot_id << std::endl;
        std::cout << "\n ok" << std::endl;
    }

    // Освобождение подсистемы MPI
    MPI_Finalize();
}
