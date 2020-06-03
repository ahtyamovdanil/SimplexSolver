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


//#define N 6
//#define M 4

//  mpiCC -o out main.cpp && mpiexec -np 8 /home/mirror/university/diploma/SimplexSolver/out

int main(int argc, char *argv[]) {

    int numprocs, myid, namelen;

    char processor_name[MPI_MAX_PROCESSOR_NAME];

    solver::init_solver(&argc, &argv, &numprocs, &myid, &namelen, processor_name);

    //int const N = 2 + 4*2;
    int M = 6;
    int N = 2 + M;
    // transposed constraint matrix
/*
    solver::matrix<float> At(N, M, {1, 1, 1, -1, 0, //x1
                                    1, -1, 2, 1, 1, //x2
                                    -1, 0, 0, 0, 0, //s1
                                    0, 1, 0, 0, 0,  //s2
                                    0, 0, 1, 0, 0,  //s3
                                    0, 0, 0, 1, 0,  //s4
                                    0, 0, 0, 0, 1,  //s5
                                    1, 0, 0, 0, 0});//t1

    std::vector<float> Cm = {-1, -1, 1, 0, 0, 0, 0, 0}; // z row
    std::vector<float> C;
    //std::vector<float> C = {-5, -4, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<float> B = {2, 1, 6, 1, 2};    // F {24, 6, 1, 2}

*/
/*
    solver::matrix<float> At(N, M, {6, 1, -1, 0,
                                    4, 2, 1, 1,
                                    1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, 1,});

*/
    solver::matrix<float> At(N, M, {1, 1, 1, -1, 0, 0,
                                    1, -1, 2, 1, 1, 1,
                                    -1, 0, 0, 0, 0, 0,
                                    0, 1, 0, 0, 0, 0,
                                    0, 0, 1, 0, 0, 0,
                                    0, 0, 0, 1, 0, 0,
                                    0, 0, 0, 0, 1, 0,
                                    0, 0, 0, 0, 0, -1});

    std::vector<float> B = {2, 1, 6, 1, 2, 0.5};
    //std::vector<float> C = {-5, -4, 0, 0, 0, 0}; // z row
    std::vector<float> C;
    //std::vector<float> B = {24, 6, 1, 2};    // F {24, 6, 1, 2}
    std::vector<float> Cm(N);


    solver::find_BFS(numprocs, myid, At, B);
    return 0; //todo
}
