//
// Created by mirror on 28.04.2020.
//

#ifndef SIMPLEXSOLVER_SOLVER_HPP
#define SIMPLEXSOLVER_SOLVER_HPP

#include <mpi.h>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>

namespace simplex{

    template<int N, int M>
    void print_matrix(std::string mess, std::array<std::array<int, M>, N> &a) {
        std::cout << mess <<std::endl;
        for(int i=0; i<N; ++i){
            for(int j=0; j<M; j++){
                std::cout << a[i][j] << ' ';
            }
            std::cout << std::endl;
        }
    };

    // initialization of MPI system
    inline void init_solver(int* argc, char** argv[], int* numprocs, int* myid, int* namelen, char* processor_name){
        MPI_Init(argc, argv);
        // Get number [numproc] of processors)
        MPI_Comm_size(MPI_COMM_WORLD, numprocs);
        // Get id [myid] of the current process
        MPI_Comm_rank(MPI_COMM_WORLD, myid);
        // Get name of the processor
        MPI_Get_processor_name(processor_name, namelen);
    }

    template<int N, int M, class T>
    void broadcast(T);
}



#endif //SIMPLEXSOLVER_SOLVER_HPP
