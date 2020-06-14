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
#include <numeric>
#include "solver.hpp"

int main(int argc, char *argv[]) {
    srand(time(nullptr));
    int N = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    std::string filename;
    bool basis = false;
    switch ( argc ){
        case 4:
            filename = argv[3];
            break;
        case 5:
            filename = argv[3];
            basis = true;
        case 3:
            basis = true;
            break;
        default:
            N = 6;
            M = 4;
    }
    int numprocs, myid, namelen;
    bool output = false;

    char processor_name[MPI_MAX_PROCESSOR_NAME];

    solver::init_solver(&argc, &argv, &numprocs, &myid, &namelen, processor_name);

    solver::matrix<float> At2;
    std::vector<float> C2(N);
    std::vector<float> B2(M);

    if(myid == 0) {
        if(argc == 3){
            At2 = solver::matrix<float>(N, M);
            std::tie(At2, B2, C2) = solver::random_problem<float>(N,M);
        }

        else{
            At2 = solver::matrix<float>(N, M);
            std::tie(At2, B2, C2) = solver::read_csv(filename, 6, 4, basis);
        }
    }

    if(!basis){
        solver::find_BFS(numprocs, myid, At2, B2);
        MPI_Finalize();
        return 0;
    }
    else{
        float optimum = 0;
        const clock_t begin_time = 0;
        bool is_rand;
        is_rand = (argc == 3);
        solver::solve(numprocs, myid, At2, C2, B2, optimum, output, is_rand);
        if(myid==0 && is_rand){
            std::cout.precision(3);
            std::cout << "time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " sec" << std::endl;
        }

        MPI_Finalize();
        return 0;
    }
}
