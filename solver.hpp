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
        matrix() : rows_(0), columns_(0) {}

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

        void add_row(std::vector<T>& row){
            if(row.size() != columns_)
                throw std::length_error("row length not equal to the number of matrix columns");
            data_.insert(data_.end(), row.begin(), row.end());
            ++ rows_;
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
        return idx;
    }


    template<typename T>
    inline std::tuple<matrix<T>, std::vector<T>, std::vector<T>>  read_csv(std::string filename){
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


    inline void solve(int& numprocs,
                      int& myid, matrix<float>& At,
                      std::vector<float>& C ,
                      std::vector<float>& B,
                      float& optimum,
                      bool output){

        int N;
        int M;

        if(myid == 0){
            N = At.rows();
            M = At.columns();
        }

        std::cout << N <<std::endl;
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::cout << "kek " << N <<std::endl;

        std::vector<int> bind(M);
        for(int i=0; i<bind.size(); ++i){
            bind[i] = i+N-M-1;
        };

        //char processor_name[MPI_MAX_PROCESSOR_NAME];
        //init_solver(argc, argv, &numprocs, &myid, &namelen, processor_name);

        // number of rows in one batch
        int batch_rows = std::floor(float(N) / numprocs);

        // batch size
        int batch_size = batch_rows * M;
        int big_batch_num = N % numprocs;
        // last batch may be smaller than the other ones

        std::vector<int> batch_sizes(numprocs);
        std::vector<int> batch_offsets(numprocs);
        std::vector<int> c_sizes(numprocs);
        std::vector<int> c_offsets(numprocs);

        for(int i=0; i<numprocs; ++i){
            if(i>=numprocs-big_batch_num){
                batch_sizes[i] = batch_size+M;
                c_sizes[i] = batch_rows+1;
                batch_offsets[i] = batch_offsets[i-1] + batch_sizes[i-1];
                c_offsets[i] = c_offsets[i-1] + c_sizes[i-1];
            }
            else{
                batch_sizes[i] = batch_size;
                c_sizes[i] = batch_rows;
                if(i==0){
                    batch_offsets[i] = 0;
                    c_offsets[i] = 0;
                }
                else{
                    batch_offsets[i] = batch_offsets[i-1] + batch_sizes[i-1];
                    c_offsets[i] = c_offsets[i-1] + c_sizes[i-1];
                }
            }
        }

        int pivot_row_idx;
        std::vector<float> pivot_column;

        // Matrix with bunch of columns for each processor
        solver::matrix<float> slave_columns(batch_rows+1, M);

        slave_columns.resize(c_sizes[myid]);

        // Vector with local pivots from all processors
        std::vector<float> local_piv_values(numprocs);
        std::vector<int> local_piv_idxs(numprocs);

        // Vector with local z row
        std::vector<float> local_c(batch_rows+1);
        local_c.resize(c_sizes[myid]);

        // indicator that optimum has been found
        bool is_done = false;

        //broadcast vector [B] to all processors
        MPI_Bcast(B.data(), M, MPI_FLOAT, 0, MPI_COMM_WORLD);

        assert(("Number of processors must be less or equal to the number of variables", N >= numprocs));
        //scatter rows of matrix [At] to each processor
        MPI_Scatterv(At.data(), batch_sizes.data(), batch_offsets.data(), MPI_FLOAT, slave_columns.data(), batch_size+M, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if(myid == 0){
            At.clear(); //todo
        }
        //scatter elements of z-row [C] to each processor
        MPI_Scatterv(C.data(), c_sizes.data(), c_offsets.data(), MPI_FLOAT, local_c.data(), batch_rows+1, MPI_FLOAT, 0, MPI_COMM_WORLD);

        std::vector<float>::iterator loc_piv_iter;
        int loc_piv_i;
        float loc_piv_val;

        std::vector<float>::iterator zvalue_iter;
        int pivot_column_proc_id;
        float global_zvalue;
        int global_zidx;

        int it = 1;

        while(!is_done){
            // find pivot column value and index in each processor
            //for(auto& item: local_c)
            //    std::cout << item << ' ';

            loc_piv_iter = std::min_element(local_c.begin(), local_c.end());
            loc_piv_i = std::distance(local_c.begin(), loc_piv_iter);
            loc_piv_val = *loc_piv_iter;
            pivot_column = slave_columns.at(loc_piv_i);

            // gather all local pivots values and indexes from all processors
            MPI_Gather(&loc_piv_val, 1, MPI_FLOAT, local_piv_values.data(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Gather(&loc_piv_i, 1, MPI_FLOAT, local_piv_idxs.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

            // get value of global pivot column and its processor id
            if (myid == 0) {
                // output for debuging
                //std::cout << "\n+-------------------------------+" << std::endl;
                //std::cout << "|  current optimum: " << optimum << "\t\t|" << std::endl;
                if(output) {
                    std::cout << "|  current optimum: " << optimum << "\t\t|" << std::endl;

                    std::cout << "|  current B: ";
                    for (int i = 0; i < B.size(); ++i)
                        std::cout << B.at(i) << '[' << bind[i] << ']' << ' ';

                    std::cout << "\t\t|\n+-------------------------------+" << std::endl;
                }
                zvalue_iter = std::min_element(local_piv_values.begin(), local_piv_values.end());
                pivot_column_proc_id = std::distance(local_piv_values.begin(), zvalue_iter);
                global_zvalue = *zvalue_iter;
                global_zidx = local_piv_idxs[pivot_column_proc_id];
                if(output){
                    std::cout << "|  current global_zvalue: " << global_zvalue << "\t\t|" << std::endl;
                    std::cout << "|  current global_zidx: " << global_zidx << "\t\t|" << std::endl;
                }
                ++it;
            }

            // broadcast processor id to all processors
            MPI_Bcast(&pivot_column_proc_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            // broadcast pivot column to all processors
            MPI_Bcast(pivot_column.data(), M, MPI_FLOAT, pivot_column_proc_id, MPI_COMM_WORLD);

            // broadcast pivot value in the z-row
            MPI_Bcast(&global_zvalue, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (global_zvalue >= 0) {
                is_done = true;
            }

            //find index of the pivot row
            if (myid == 0) {
                pivot_row_idx = solver::find_piv_row(B, pivot_column);
                if(pivot_row_idx == -1)
                    break;
                int it = std::accumulate(c_sizes.begin(), c_sizes.begin()+pivot_column_proc_id, 0);
                if(output)
                    std::cout << "|  current pivot_row_idx: " << pivot_row_idx << "\t\t|" << std::endl;
                bind[pivot_row_idx] = it + global_zidx;
            }

            MPI_Bcast(&pivot_row_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // pivot value = element on the cross of the pivot row and pivot column
            float pivot_elem = pivot_column[pivot_row_idx];
            // update pivot row in each processor
            for (int i=0; i<slave_columns.rows(); ++i) {
                slave_columns.at(i,pivot_row_idx) /= pivot_elem;
            }
            // update [B] in pivot row
            B[pivot_row_idx] /= pivot_elem;

            // update global optimum value
            if (myid == 0) {
                optimum -= global_zvalue * B[pivot_row_idx];
            }

            // update z-row
            for (int j = 0; j < local_c.size(); ++j) {
                auto abc = slave_columns.at(j,pivot_row_idx);
                local_c[j] -= global_zvalue * slave_columns.at(j,pivot_row_idx);
            }

            //update values in the rest of the table
            for (int i=0; i<slave_columns.rows(); ++i) {
                for (int j = 0; j < M; ++j) {
                    if (j != pivot_row_idx) {
                        slave_columns.at(i,j) -= pivot_column[j] * slave_columns.at(i,pivot_row_idx);
                    }
                }
            }

            // update results vector [B]
            for (int j = 0; j < M; ++j) {
                if (j != pivot_row_idx) {
                    B[j] -= pivot_column[j] * B[pivot_row_idx];
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

        // gather vector C from all processors
        MPI_Gatherv(local_c.data(), batch_rows,  MPI_FLOAT, C.data(),c_sizes.data(), c_offsets.data(), MPI_FLOAT, 0, MPI_COMM_WORLD);

        // gather all columns
        //MPI_Gather(slave_columns.data(), batch_size, MPI_FLOAT, At.data(), batch_size, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        // print finale result
        if(myid==0 && output){
            std::cout << "----------------------" <<std::endl;
            std::cout << "\nB:" << std::endl;

            for(int i=0; i<B.size(); ++i)
                std::cout << B[i] << ' ';
            //std::cout << "\nC:" << std::endl;
            if(is_done)
                std::cout << "\noptimum value: " << optimum << std::endl;
            else
                std::cout << "\nproblem is unbounded: "<< std::endl;
            std::cout << "number of iterations: " << it << std::endl;
            for(int i=0; i < bind.size(); ++i){
                std::cout << bind[i] << ' ';
            }
            std::cout << "\ndone" << std::endl;
        }
    }

    inline void find_BFS(int& numprocs, int& myid, matrix<float>& At, std::vector<float>& B) {
        int N = At.rows();
        int M = At.columns();

        std::vector<int> bind(M);
        std::vector<float> C(N);

        for (int i = 0; i < bind.size(); ++i) {
            bind[i] = i + N - M - 1;
        };
        float optimum = 0;

        //add rows to the matrix At
        for(int i=0; i<M; ++i){
            if(std::signbit(At.at(i+N-M, i)) != std::signbit(B[i])){
                C.push_back(0);
                for(int j=0; j<N; ++j)
                    C[j] -= At.at(j, i);
                optimum -= B[i];
                bind[i] = N+i;
                std::vector<float> row(M);
                row[i] = 1;
                At.add_row(row);
            }
        }

        solve(numprocs, myid, At, C, B, optimum, false);
    }
}

#endif //SIMPLEXSOLVER_SOLVER_HPP
