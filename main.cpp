#define N 6
#define M 4

#include <mpi.h>
#include <iostream>
#include <vector>
#include <array>
#include <cassert>


//  mpiCC -o out main.cpp && mpiexec -np 8 /home/mirror/university/diploma/SimplexSolver/out
void print_results(std::string mess, std::array<std::array<int, M>, N> a){
    std::cout << mess <<std::endl;
    for(int i=0; i<N; ++i){
        for(int j=0; j<M; j++){
            std::cout << a[i][j] << ' ';
        }
        std::cout << std::endl;
    }
};

int main(int argc, char *argv[])
{
    int numprocs, myid, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    // Инициализация подсистемы MPI
    MPI_Init(&argc, &argv);
    // Получить общее число процессов numprocs в рамках задачи)
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // Получить номер myid текущего процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);

    //std::cout << "number of processors: " << numprocs << std::endl;

    using Matrix = std::array<std::array<int, N>, M>;
    using MatrixT = std::array<std::array<int, M>, N>;

    int batch_rows = N / numprocs;
    int  batch_size  =  batch_rows * M;
    // constraint matrix
    Matrix A = {{{6,4,1,0,0,0},
                 {1,2,0,1,0,0},
                 {-1,1,0,0,1,0},
                 {0,1,0,0,0,1}}};

    // transposed constraint matrix
    MatrixT At = {{{6,1,-1,0},
                   {4,2,1,1},
                   {1,0,0,0},
                   {0,1,0,0},
                   {0,0,1,0},
                   {0,0,0,1}}};

    MatrixT At_test;
    std::array<int, N> C = {-5,-4,0,0,0,0}; // z row
    std::array<int, M> B = {24,6,1,2};      // F(x)

    std::vector<std::array<int, M>> slave_columns(batch_rows);

    //broadcast vector C to all processes
    MPI_Bcast(C.data(), N, MPI_INT, 0, MPI_COMM_WORLD);
    //broadcast vector B to all processes
    MPI_Bcast(B.data(), M, MPI_INT, 0, MPI_COMM_WORLD);

    //scatter rows of matrix At to different processes
    //ToDo
    assert(("Number of processors must be a multiple of the number of variables", N % numprocs == 0));

    MPI_Scatter(At.data(), batch_size, MPI_INT, slave_columns.at(0).data(), batch_size, MPI_INT, 0 , MPI_COMM_WORLD);

    for(auto &col: slave_columns){
        for(auto &item: col) {
            item = myid;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(slave_columns.at(0).data(), batch_size, MPI_INT, At_test.data(), batch_size, MPI_INT, 0, MPI_COMM_WORLD);

    if(myid==0)
        print_results("results", At_test);

    // Освобождение подсистемы MPI
    MPI_Finalize();
}
/*array([[ 54,  37,  47,  57],
       [130,  93, 119, 145],
       [ 44,  41,  56,  71],
       [111,  79, 101, 123]])
*/
