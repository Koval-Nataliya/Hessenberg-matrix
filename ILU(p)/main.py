import scipy.sparse as sp
from mpi4py import MPI
import numpy as np


class SparseMatrix:
    def __init__(self, rows, columns):
        self.values = []
        self.row_indices = [0] * (rows + 1)
        self.column_indices = []
        self.rows_count = rows
        self.columns_count = columns

    def __getitem__(self, index):
        row, column = index
        n1 = self.row_indices[row + 1]
        for k in range(self.row_indices[row], n1):
            if self.column_indices[k] == column:
                return self.values[k]
        return 0

    def __setitem__(self, index, value):
        row, column = index
        n1 = self.row_indices[row + 1]
        for k in range(self.row_indices[row], n1):
            if self.column_indices[k] == column:
                self.values[k] = value
                return
        k = self.row_indices[row + 1]
        self.values.insert(k, value)
        self.column_indices.insert(k, column)
        for i in range(row + 1, len(self.row_indices)):
            self.row_indices[i] += 1

    def __call__(self, row, column):
        return self.__getitem__((row, column))

    def __str__(self):
        output = ''
        for i in range(self.rows_count):
            for j in range(self.columns_count):
                output += str(self[i, j]) + ' '
            output += '\n'
        return output

class SparseVector:

    def __init__(self, size=0):
        self.size = size
        self.values = []
        self.indices = []

    def get_size(self):
        return len(self.indices)

    def __getitem__(self, index):
        if index >= self.size:
            raise Exception("Index out of range")

        ind = self.indices.index(index) if index in self.indices else -1

        if ind != -1:
            return self.values[ind]
        else:
            return 0

    def __setitem__(self, index: int, value: float):
        if index >= self.size:
            raise Exception("Index out of range")

        ind = self.indices.index(index) if index in self.indices else -1
        if ind == -1:
            self.indices.append(index)
            self.values.append(value)
        else:
            self.values[ind] = value

    def __eq__(self, other):
        self.values = other.values
        self.indices = other.indices
        self.size = other.size
        return self


def ILUp(matrix, n, p):
    lev = [0]*(n*n)
    cw = MPI.COMM_WORLD
    num_tasks = cw.Get_size()
    task_id = cw.Get_rank()
    start_time = MPI.Wtime()

    mat = SparseMatrix(n, n)
    number_of_proc = 4

    block_size = n // number_of_proc

    chunk_matrix = SparseMatrix(n, n)
    if task_id == 0:
        mat.values = matrix.values.copy()
        mat.row_indices = matrix.row_indices.copy()
        mat.column_indices = matrix.column_indices.copy()
        mat.rows_count = matrix.rows_count
        mat.columns_count = matrix.columns_count
        for i in range(n):
            for j in range(n):
                if mat[i, j] != 0 or (i == j):
                    lev[i*n+j] = 0
                else:
                    lev[i*n + j] = p+1
        for i in range(1, number_of_proc):
            first_r = i * block_size
            if i == number_of_proc - 1:
                last_r = n
            else:
                last_r = (i + 1) * block_size

            rowPtr_tmp = np.array([],dtype=int)
            values_tmp = np.array([], dtype=float)
            colInd_tmp = np.array([], dtype=int)

            for j in range(first_r, last_r + 1):
                rowPtr_tmp = np.append(rowPtr_tmp, mat.row_indices[j] - mat.row_indices[first_r])

            start = mat.row_indices[first_r]
            end = mat.row_indices[last_r]

            for j in range(start, end):
                values_tmp = np.append(values_tmp, mat.values[j])
                colInd_tmp = np.append(colInd_tmp, mat.column_indices[j])

            s1 = last_r - first_r + 1
            s2 = end - start
            s3 = end - start

            size = np.array([s1, s2, s3], dtype=int)
            cw.Send(size, dest=i, tag=0)
            rowPtr_tmp =np.array(rowPtr_tmp, dtype=int)
            colInd_tmp =colInd_tmp.astype(int)
            values_tmp =values_tmp.astype(float)

            cw.Send(rowPtr_tmp, dest=i, tag=1)
            cw.Send(colInd_tmp, dest=i, tag=2)
            cw.Send(values_tmp, dest=i, tag=3)
            cw.Send(lev, dest=i, tag=4)


        for j in range(block_size + 1):
            rowPtr_tmp = np.append(rowPtr_tmp, mat.row_indices[j])

        start = mat.row_indices[0]
        end = mat.row_indices[block_size]

        for j in range(start, end):
            values_tmp = np.append(values_tmp, mat.values[j])
            colInd_tmp = np.append(colInd_tmp, mat.column_indices[j])

        chunk_matrix.column_indices = colInd_tmp.copy()
        chunk_matrix.row_indices = rowPtr_tmp.copy()
        chunk_matrix.values = values_tmp.copy()
        chunk_matrix.rows_count = block_size
        chunk_matrix.columns_count = n

    if task_id < number_of_proc and task_id != 0:
        size = np.empty(3 + 1, dtype=int)
        cw.Recv(size, source=0, tag=0)

        rowPtr_tmp = np.empty(size[0] + 1 , dtype=int)
        cw.Recv(rowPtr_tmp, source=0,tag=1)

        colInd_tmp = np.empty(size[1] + 1, dtype=int)
        cw.Recv(colInd_tmp, source=0, tag=2)

        values_tmp = np.empty(size[2] + 1, dtype=float)
        cw.Recv(values_tmp, source=0, tag=3)
        cw.Recv(lev, source=0, tag=4)

        chunk_matrix.column_indices = colInd_tmp.copy()
        chunk_matrix.row_indices = rowPtr_tmp.copy()
        chunk_matrix.values = values_tmp.copy()
        chunk_matrix.rows_count = size[0] - 1
        chunk_matrix.columns_count = n

    tmp = np.zeros(n)
    result = SparseMatrix(n, n)
    count = 0

    for i in range(n):
        if i // block_size == task_id or (i // block_size > number_of_proc - 1 and num_tasks == number_of_proc - 1):
            local_row = i % block_size
            for k in range(i):
                if k >= count:
                    r = np.empty(n + 1, dtype=float)
                    cw.Recv(r, source=MPI.ANY_SOURCE, tag=0)

                    for j in range(n):
                        result[count, j] = r[j]
                    count += 1

                if result[k, k] != 0 and tmp[k] != 0:
                    lev[local_row*n + k] = lev[local_row*n + local_row] + lev[k*n + local_row] + 1
                    if lev[n * local_row + k] > p:
                        tmp[k] = 0
                    else:
                        tmp2 = tmp[k]/result[k,k]
                        for j in range(result.row_indices[k] + 1, result.row_indices[k + 1]):
                            if result.column_indices[j] >= k + 1:
                                tmp[result.column_indices[j]] -= tmp2 * result.values[j]

            for j in range(task_id + 1, number_of_proc):
                req  = cw.Isend(tmp, dest=j, tag=0)
                req.Wait()

            for j in range(n):
                result[i, j] = tmp[j]
            count += 1
            tmp.fill(0)

    if task_id == number_of_proc - 1:
        end_time = MPI.Wtime()
        print("Program execution time: ", (end_time - start_time), " seconds")
    return result


n = 200

np.random.seed(42)
A = np.random.choice(a=[0, 1, 2, 3, 4, 5], p=[0.5, 0.1, 0.1, 0.1, 0.1, 0.1], size=(n, n))
A = sp.csr_matrix(A)

p = 2

m = SparseMatrix(n, n)
m.values = list(A.data)
m.row_indices = list(A.indptr)
m.column_indices = list(A.indices)
m.rows_count = n
m.columns_count = n

ILUp(m, n, p)

