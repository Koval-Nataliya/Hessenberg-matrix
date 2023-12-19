#include <iostream>
#include <vector>
#include <algorithm>    // std::find
#include <utility>
#include <map>
#include "mpi.h"
int p = 3;
struct SparseVec {
    int nz;
    std::vector<double> ind;
    std::vector<double> val;
    std::vector<int> levels;
};
struct SparseVec axy(SparseVec v1, SparseVec v2, int i) {
    //std::cout<<"in func"<<'\n';

    struct SparseVec tmp;
    double a1, a2, coef;
    auto f2 =std::find(v2.ind.begin(), v2.ind.end(), i);
    //если нет элемента, который надо занулить
    if (f2 == v2.ind.end()){
        return v2;
    } else {
        //find coef
        auto f1 =std::find(v1.ind.begin(), v1.ind.end(), i);
        int t = f1- v1.ind.begin();
        a1 = v1.val[t];
        t = f2 - v2.ind.begin();
        a2 = v2.val[t];
        coef = a2/a1;
        int i1 = 0, i2 = 0;
        while (i1 < v1.nz && i2 < v2.nz) {
            //std::cout<<i1<<' '<<i2<<'\n';
            if(v1.ind[i1] == v2.ind[i2]) {
                tmp.ind.push_back(i1);
                tmp.val.push_back(v2.val[i2] - coef*v1.val[i1]);
                tmp.levels.push_back(std::max(v1.levels[i1], v2.levels[i2]));
                i1++; i2++;
            } else if (v1.ind[i1] < v2.ind[i2]) {
                if(v1.levels[i1] < p) {
                    tmp.ind.push_back(i1);
                    tmp.val.push_back(v1.val[i1]);
                    tmp.levels.push_back(v1.levels[i1]+1);
                    i1++;
                }
            } else {
                if(v2.levels[i2] < p) {
                    tmp.ind.push_back(i2);
                    tmp.val.push_back(v2.val[i2]);
                    tmp.levels.push_back(v2.levels[i2]+1);
                    i2++;
                }
            }
            i1++; i2++;
        }
        std::sort(tmp.ind.begin(), tmp.ind.end());
        tmp.nz = tmp.ind.size();
    }
    return tmp;
}

int main() {
    //Initialize MPI
    MPI_Init(NULL, NULL);
    // Get total number of tasks

    int number_of_str = 10;
    struct SparseVec vec[number_of_str];
    int n = 0;
    //filling matrix;
    for (int i = 0; i < number_of_str; i++) {
        n = rand() % 5;
        vec[i].nz = n;
        //чтобы был диагональный элемент
        vec[i].ind.push_back(i);
        vec[i].val.push_back(rand() % 100);
        vec[i].levels.push_back(1);
        //остальные элементы
        for (int j = 0; j < n - 1; j++) {
            vec[i].ind.push_back(vec[i].nz - int((j - i) / 2) + rand() % 4);
            vec[i].val.push_back(rand() % 100);
            vec[i].levels.push_back(1);
        }
        std::sort(vec[i].ind.begin(), vec[i].ind.end());
    }

    int num_tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    // Get the task ID
    int task_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    double start_time, end_time;
    start_time = MPI_Wtime();

    int number_of_proc = 4;

    // matrix block size for each process
    int block_size = number_of_str / number_of_proc;

    MPI_Datatype my_struct_type;
    int blocklengths[4] = {1, number_of_str, number_of_str, number_of_str};
    MPI_Datatype types[4] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Aint displacements[4] = {offsetof(my_struct, nz), offsetof(my_struct, ind), offsetof(my_struct, val),
                                 offsetof(my_struct, levels)};
    MPI_Type_create_struct(4, blocklengths, displacements, types, &my_struct_type);
    MPI_Type_commit(&my_struct_type);
    my_struct vec;
    struct SparseVec mat_n[block_size];
    if (task_id == 0) {
        for (int i = 0; i < number_of_proc; i++) {
            const int start_row = i * block_size;
            int end_row;
            if (i == number_of_proc - 1) {
                end_row = n;
            } else {
                end_row = (i + 1) * block_size;
            }
            for (int j = start_row; j < end_row + 1; j++) {
                mat_n[j - start_row] = vec[j];
            }
            for (int k = 0; k < end_row; k++) {
                for (int j = i + 1; j < end_row; j++) {
                    mat_n[j] = axy(mat_n[k], mat_n[j], k);
                    vec[j] = axy(vec[k], vec[j], k);
                }
            }
            MPI_Send(&mat_n, block_size, my_struct_type, 1, 0, MPI_COMM_WORLD);
            MPI_Send(start_row, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
            MPI_Send(end_row, 1, MPI_INT, 1, 2, MPI_COMM_WORLD);
        }
    }
    if (task_id < number_of_proc && task_id != 0) {
        MPI_Recv(&mat_n, block_size, my_struct_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(start_row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(end_row, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        for (int i = end_row; i < number_of_str; i++) {
            for (int j = 0; j < end_row; j++) {
                vec[i] = axy(mat_n[j], vec[i], i);
            }
        }
    }

    return 0;
}
