#define N 4

#include <mpi.h>
#include <iostream>
#include <vector>
#include <array>

void print_results(std::string mess, std::array<std::array<int, N>, N> a){
    std::cout << mess <<std::endl;
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; j++){
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

    using Matrix = std::array<std::array<int, N>, N>;
    Matrix a = {{{1,2,3,4},{5,6,7,8},{9,1,2,3},{4,5,6,7}}};
    Matrix b = {{{1,2,3,4},{5,6,7,8},{9,1,2,3},{4,5,6,7}}};
    Matrix c;
    std::array<int, N> row_vec{}, cc{};
    bool done = false;
    int sum;
    int vec_count = numprocs;
    int n_batches = N / numprocs;
    int curr_row=0;
    //broadcast second matrix to all processes
    MPI_Bcast(b.data(), N*N, MPI_INT, 0, MPI_COMM_WORLD);

    //scatter rows of first matrix to different processes
    for(;curr_row < n_batches*vec_count; curr_row+=vec_count){
        MPI_Scatter(&a[curr_row], N, MPI_INT, row_vec.data(), N, MPI_INT, 0 , MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        //perform vector multiplication by all processes
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                sum = sum + row_vec[j] * b[j][i];
            }
            cc[i] = sum;
            sum = 0;
        }
        MPI_Gather(cc.data(), N, MPI_INT, &c[curr_row], N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(N % numprocs != 0){
        std::array<int, 1000> scounts = {0};
        std::array<int, 1000> displs = {0};
        for(int k=0; k < N % numprocs; ++k){
            scounts[k] = N;
            displs[k] = k*N;
        }

        //perform vector multiplication by all processes
        if(scounts[myid]!=0){
            MPI_Scatterv(&a[curr_row], scounts.data(), displs.data(), MPI_INT, row_vec.data(), N, MPI_INT, 0 , MPI_COMM_WORLD);
            for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                    sum = sum + row_vec[j] * b[j][i];
                }
                cc[i] = sum;
                sum = 0;
            }
        }
        MPI_Gatherv(cc.data(), N, MPI_INT, &c[curr_row], scounts.data(), displs.data() ,MPI_INT, 0, MPI_COMM_WORLD);
    }

    if(myid==0)
        print_results("results", c);

    // Освобождение подсистемы MPI
    MPI_Finalize();
}
/*array([[ 54,  37,  47,  57],
       [130,  93, 119, 145],
       [ 44,  41,  56,  71],
       [111,  79, 101, 123]])
*/
