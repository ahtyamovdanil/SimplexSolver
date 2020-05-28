//
// Created by Ahtyamov Danil on 28.04.2020.
// https://github.com/ahtyamovdanil
//
#ifndef SIMPLEXSOLVER_SOLVER_HPP
#define SIMPLEXSOLVER_SOLVER_HPP

#include <mpi.h>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace simplex {

    template<class T>
    void print_matrix(std::string mess, T &a) {
        std::cout << mess << std::endl;
        for (auto & row: a) {
            for (auto & item: row) {
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }
    };

    // initialization of MPI system
    inline void init_solver(int *argc, char **argv[], int *numprocs, int *myid, int *namelen, char *processor_name) {
        MPI_Init(argc, argv);
        // Get number [numproc] of processors)
        MPI_Comm_size(MPI_COMM_WORLD, numprocs);
        // Get id [myid] of the current process
        MPI_Comm_rank(MPI_COMM_WORLD, myid);
        // Get name of the processor
        MPI_Get_processor_name(processor_name, namelen);
    }

    template<class T>
    int find_piv_row(T &arr1, T &arr2) {
        if (arr1.size() != arr2.size())
            throw std::range_error("size of arrays must be same");
        int idx, result, val = -1;
        for (int i = 0; i < arr1.size(); i++) {
            if (arr2[i] != 0) {
                val = arr1[i] / arr2[i];
                if (val > 0 && (val < result || result == -1)) {
                    std::cout << "res: " << val << std::endl;
                    result = val;
                    idx = i;
                }
            }
        }
        std::cout << "result " << result << std::endl;
        return idx;
    }

}

#endif //SIMPLEXSOLVER_SOLVER_HPP
