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
#include <tuple>
#include <ctime>

namespace solver {

    template <class T>
    class matrix {
        std::vector<T> data_;
        int rows_;
        int columns_;
    public:
        matrix(int rows, int columns) :rows_(rows), columns_(columns) {
            data_.reserve(rows*columns);
        }
        matrix(int rows, int columns, std::vector<T> data) : rows_(rows), columns_(columns), data_(data) {}

        T &operator()(int row, int column) { return data_[row*columns_+column]; }
        T &operator()(int index) { return data_[index]; }

        int rows(){return rows_;}
        int columns(){return columns_;}
        T* data(){ return data_.data();}

        void resize(int rows){
            rows_ = rows;
        };

        void clear(){
            rows_ = 0;
            columns_ = 0;
            data_.clear();
        }

        T& at(int row, int column){
            return data_[row*columns_+column];
        }
        std::vector<T> at(int row){
            return std::vector<T>(data_.begin() + row*columns_, data_.begin() + (row+1)*columns_);
        }


        //~matrix(){ std::cout<<"destructor\n"; data_.~vector(); }
    };

    template<typename T>
    void print_matrix(std::string mess, T &a) {
        std::cout << mess << std::endl;
        for (int i=0; i<a.size(); ++i) {
            for (int j=0; j<a.at(i).size(); ++j) {
                std::cout << a[i][j] << " ";
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

    inline void solve(matrix<float>& At, std::vector<float>& C , std::vector<float>& B){

    }

    template<typename T>
    int find_piv_row(T &B, T &piv_col) {
        if (B.size() != piv_col.size()){}
            //throw std::range_error("size of arrays must be same");
        int idx=-1;
        float result=-1, val = -1;
        for (int i = 0; i < B.size(); ++i) {
            if (piv_col[i] != 0) {
                val = B[i] / piv_col[i];
                if (val > 0 && (val < result || result == -1)) {
                    //std::cout << "res: " << val << std::endl;
                    result = val;
                    idx = i;
                }
            }
        }
        if (idx == -1)
            throw std::out_of_range("kek");
        //std::cout << "result " << result << std::endl;
        return idx;
    }

    inline void read_csv(std::string filename){
        // todo
    }

    template<typename T>
    inline std::tuple<matrix<T>, std::vector<T>, std::vector<T>> random_problem(int N, int M) {
        srand(time(NULL));
        matrix<T> At(N, M);
        std::vector<T> B(M);
        std::vector<T> C(N);
        //std::vector<T> expected(N);
        T min = -1;
        T max = 1;
        for(int i=0; i<N; ++i){
            for(int j=0; j<M; ++j){
                //At.at(i,j) = min + (((float) rand()) / (float) RAND_MAX) * (max - min);
                if(i >= N-M){
                    if(i-N+M == j)
                        At.at(i,j) = 1;
                    else
                        At.at(i,j) = 0;
                }
                else
                    At.at(i,j) = rand() % 21 - 4;
            }
        }

        for(int i=0; i<N; ++i){
            if(i<N-M)
                C[i] = -1*(rand() % 10 + 1);
                //C[i] = min + (((float) rand()) / (float) RAND_MAX) * (max - min);
            else
                C[i] = 0;
        }

        for(int i=0; i<M; ++i)
            B[i] = rand() % 10 + 1;

/*
        T eta = 0.1;
        for(int i=0; i<M; ++i){
            T bi = 0;
            for(int j=0; j<N; ++j){
                bi += At.at(j,i)*expected[j];
            }
            B[i] = bi + eta;
        }
        for(int i=0; i<N; ++i)
            std::cout << expected[i] << ' ';
        std::cout<<std::endl<<std::endl;
        */
/*
        for(int i=0; i<M; ++i){
            for(int j=0; j<N; ++j){
                std::cout << At.at(j, i) << ' ';
            }
            std::cout << std::endl;
        }
        for(int j=0; j<M; ++j){
            std::cout << B.at(j) << ' ';
        }
        std::cout << std::endl;

        for(int j=0; j<N; ++j){
            std::cout << C.at(j) << ' ';
        }
        std::cout << std::endl;
*/
        return std::make_tuple(std::move(At),std::move(B),std::move(C));
    }

}

#endif //SIMPLEXSOLVER_SOLVER_HPP
