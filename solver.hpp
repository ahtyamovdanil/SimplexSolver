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

    template<typename T>
    int find_piv_row(T &arr1, T &arr2) {
        if (arr1.size() != arr2.size()){}
            //throw std::range_error("size of arrays must be same");
        int idx, result, val = -1;
        for (int i = 0; i < arr1.size(); i++) {
            if (arr2[i] != 0) {
                val = arr1[i] / arr2[i];
                if (val > 0 && (val < result || result == -1)) {
                    //std::cout << "res: " << val << std::endl;
                    result = val;
                    idx = i;
                }
            }
        }
        //std::cout << "result " << result << std::endl;
        return idx;
    }

    inline void read_csv(std::string filename){

    }

}

#endif //SIMPLEXSOLVER_SOLVER_HPP
