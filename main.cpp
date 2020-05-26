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

    using Matrix = std::array<std::array<int, N>, M>;
    using MatrixT = std::array<std::array<int, M>, N>;

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
    std::array<int, N> C = {-5,-4,0,0,0}; // z row
    std::array<int, M> B = {24,6,1,2};      // F(x)

    std::vector<std::array<int, M>> slave_columns(batch_rows);

    //broadcast vector C to all processes
    MPI_Bcast(C.data(), N, MPI_INT, 0, MPI_COMM_WORLD);
    //broadcast vector B to all processes
    MPI_Bcast(B.data(), M, MPI_INT, 0, MPI_COMM_WORLD);

    //scatter rows of matrix At to different processes
    assert(("Number of processors must be less or equal to the number of variables", N >= numprocs));

    MPI_Scatter(At.data(), batch_size, MPI_INT, slave_columns.at(0).data(), batch_size, MPI_INT, 0 , MPI_COMM_WORLD);

    // last process may have smaller batch then other ones
    if(myid == numprocs - 1){
        slave_columns.resize(rows_in_last);
    }

    for(auto &col: slave_columns){
        for(auto &item: col) {
            item = myid;
            //std::cout << myid << ' ' << item << std::endl;
        }
    }

    MPI_Gather(slave_columns.at(0).data(), batch_size, MPI_INT, At_test.data(), batch_size, MPI_INT, 0, MPI_COMM_WORLD);

    if(myid==0){
        simplex::print_matrix<N, M>("results", At_test);
        std::cout << "ok" << std::endl;
    }

    // Освобождение подсистемы MPI
    MPI_Finalize();
}
