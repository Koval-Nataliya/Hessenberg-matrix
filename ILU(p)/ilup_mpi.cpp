#include <vector>
#include <cmath>
#include <iostream>
#include <chrono>
#include "mpi.h"
#include <fstream>

template<class T>
class SparseMatrix
{
public:
    std::vector<T> val;
    std::vector<int> row_ind;
    std::vector<int> column_ind;
    int rows_count, columns_count;

    SparseMatrix() : rows_count(0), columns_count(0) {}
    ~SparseMatrix() = default;
    SparseMatrix(int rows, int columns)
        : val(),
          row_ind(rows + 1),
          column_ind(),
          rows_count(rows),
          columns_count(columns) {}


    double& operator()(int row, int column)

    {
        int n1 = row_ind[row + 1];

        for (int k = row_ind[row]; k < n1; k++)
            if (column_ind[k] == column) return val[k];
        int k = row_ind[row + 1];
        val.insert(val.begin() + k, 0);
        column_ind.insert(column_ind.begin() + k, column);

        for (int i = row + 1, n2 = row_ind.size(); i < n2; i++)
            row_ind[i]++;
        return val[k];
    }

    double operator()(int row, int column) const
    {
        int n = row_ind[row + 1];
        for (int k = row_ind[row]; k < n; ++k)
            if (column_ind[k] == column)
                return val[k];
        return 0;
    }
};

SparseMatrix<double>ILUp(SparseMatrix<double> matrix, int n, double p)
{
    std::vector <int> lev(n*n);

    int num_tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    int task_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    double start_time, end_time;
    start_time = MPI_Wtime();
    SparseMatrix<double> mat(n, n);
    int number_of_proc = 4;
    int block_size = n / number_of_proc;

    SparseMatrix<double> chunk_matrix;
    SparseMatrix<double> res;

    std::vector<int> size;
    std::vector<int> rowPtr_tmp;
    std::vector<int> colInd_tmp;
    std::vector<double> values_tmp;

    if (task_id == 0)
    {
        mat.row_ind = matrix.row_ind;
        mat.column_ind = matrix.column_ind;
        mat.val = matrix.val;
        mat.rows_count = matrix.rows_count;
        mat.columns_count = matrix.columns_count;

        for(int i =0; i < n; i++) {
            for (int j = 0;j < n; j++) {
                if (mat(i, j) != 0 or (i ==j)) {
                    lev[i*n + j] = 0;
                } else{
                    lev[i*n + j] = p+1;
                }

            }
        }

        for (int i = 1; i < number_of_proc; i++)
        {
            const int first_r = i * block_size;
            int last_r;
            if (i == number_of_proc - 1)
            {
                last_r = n;
            }
            else {
                last_r = (i + 1) * block_size;
            }

            for (int j = first_r; j < last_r + 1; j++)
            {
                rowPtr_tmp.push_back(mat.row_ind[j] - mat.row_ind[first_r]);
            }

            auto start = mat.row_ind[first_r];
            auto end = mat.row_ind[last_r];

            for (int j = start; j < end; j++)
            {
                values_tmp.push_back(mat.val[j]);
                colInd_tmp.push_back(mat.column_ind[j]);
            }

            int s1 = last_r - first_r + 1;
            int s2 = end - start;
            int s3 = end - start;

            size = { s1, s2, s3 };
            MPI_Send(size.data(), 3, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(lev, n*n, MPI_INT, i, 1, MPI_COMM_WORLD);

            MPI_Send(rowPtr_tmp.data(), s1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(colInd_tmp.data(), s2, MPI_INT, i, 2, MPI_COMM_WORLD);
            MPI_Send(values_tmp.data(), s3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
        }

        for (int j = 0; j < block_size + 1; j++)
        {
            rowPtr_tmp.push_back(mat.row_ind[j]);
        }
        auto start = mat.row_ind[0];
        auto end = mat.row_ind[block_size];

        for (int j = start; j < end; j++)
        {
            values_tmp.push_back(mat.val[j]);
            colInd_tmp.push_back(mat.column_ind[j]);
        }

        chunk_matrix.column_ind = colInd_tmp;
        chunk_matrix.row_ind = rowPtr_tmp;
        chunk_matrix.val = values_tmp;
        chunk_matrix.rows_count = block_size;
        chunk_matrix.columns_count = n;

    }

    if (task_id < number_of_proc && task_id != 0)
    {

        size.resize(3);
        MPI_Recv(size.data(), 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        rowPtr_tmp.resize(size[0]);
        MPI_Recv(rowPtr_tmp.data(), size[0], MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(lev, n*n, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        colInd_tmp.resize(size[1]);
        MPI_Recv(colInd_tmp.data(), size[1], MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        values_tmp.resize(size[2]);
        MPI_Recv(values_tmp.data(), size[2], MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        chunk_matrix.column_ind = colInd_tmp;
        chunk_matrix.row_ind = rowPtr_tmp;
        chunk_matrix.val = values_tmp;
        chunk_matrix.rows_count = size[0] - 1;
        chunk_matrix.columns_count = n;
    }
    std::vector<double> tmp(n, 0);
    SparseMatrix<double> result(n, n);

    int count = 0;

    for (int i = 0; i < n; i++)
    {
        if (i / block_size == task_id || i / block_size > number_of_proc - 1 && task_id == number_of_proc - 1)
        {
            int local_row = i % block_size;

            for (int k = 0; k < i; k++)
            {
                if (result(k, k) != 0 && tmp[k] != 0)
                {
                    lev[local_row*n + k] = lev[local_row*n + local_row] + lev[k*n + local_row] + 1;
                    if (lev[n*local_row + k] > p)
                    {
                        tmp[k] = 0;
                    }
                    else {
                        const double tmp2 = (tmp[k] / result(k, k));

                        for (int j = result.row_ind[k] + 1; j < result.row_ind[k + 1]; j++)
                        {
                            if (result.column_ind[j] >= k + 1)
                            {
                                tmp[result.column_ind[j]] -= tmp2 * result.val[j];
                            }
                        }
                    }
                }
            }
            std::vector<MPI_Request>rec(number_of_proc - task_id);
            for (int t = task_id + 1; t < number_of_proc; t++)
            {
                MPI_Isend(tmp.data(), n, MPI_DOUBLE, t, 0, MPI_COMM_WORLD, &rec[t - task_id - 1]);
                MPI_Wait(&rec[t - task_id - 1], MPI_STATUS_IGNORE);
            }
            for (int j = 0; j < n; j++)
            {
                result(i, j) = tmp[j];
            }
            count++;
            tmp.clear();
            tmp.resize(n, 0);
        }
    }
    if (task_id == number_of_proc - 1)
    {
        end_time = MPI_Wtime();
        std::cout << "Program execution time: " << (end_time - start_time) << " seconds" << std::endl;
    }
    return result;
}


int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    int n = 10;
    double p;
    std::cin>>p;
    SparseMatrix<double> matr(n, n);

    std::vector<int> rowPtr;
    std::vector<int> colInd;
    std::vector<double> values;

    matr.val = {2, 4, 2, 4, 1, 1, 1, 4, 1, 8, 4, 1, 2, 2, 5, 1, 4, 1, 1, 4, 6, 4, 5, 7, 1, 3};
    matr.column_ind = { 4, 0, 1, 5, 0, 6, 7, 3, 7, 0, 3, 4, 7, 0, 1, 3, 4, 7, 0, 1, 3, 0, 3, 4, 6, 7 };
    matr.row_ind = {0, 1, 4, 7, 9, 13, 18, 21, 26};
    matr.rows_count = n;
    matr.columns_count = n;

    ILUp(matr, n, p);
    MPI_Finalize();

    return 0;
}
