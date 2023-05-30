#include <vector>
#include <iostream>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <exception>
#include <omp.h>

int num_threads = 6;

template<class T>
class SparseVector {
private:
    size_t size;
public:
    std::vector<T> values;       // values of non-zero elements
    std::vector<size_t> indices; // indexes of non-zero elements

    SparseVector() :size(0) {}
    ~SparseVector() = default;

    SparseVector(size_t size): values(0), indices(0), size(size) {}

    size_t get_size() const { return indices.size(); }

    double& operator[](size_t index)
    {
        if (index >= size) { throw std::exception(); }

        auto ind = std::find(indices.begin(), indices.end(), index);
        if (ind != indices.end())
        {
            return values[ind - indices.begin()];
        }
        else
        {
            indices.push_back(index);
            values.push_back(0);
            return values[values.size() - 1];
        }
    }

    double operator[](size_t index) const
    {
        if (index >= size) { throw std::exception(); }

        auto ind = std::find(indices.begin(), indices.end(), index);
        if (ind == indices.end())
        {
            return 0;
        }
        else
        {
            return values[ind - indices.begin()];
        }
    }

    // Overloading the equals operator to assign values
    SparseVector<T>& operator=(const SparseVector<T>& other)
    {
        values = other.values;
        indices = other.indices;
        return *this;
    }

    void zero_reset()
    {
        values.clear();
        indices.clear();
    }

};


template<class T>
class SparseMatrix
{
public:
    std::vector<T> values;                  
    std::vector<size_t> row_indices;        
    std::vector<size_t> column_indices;     
    size_t rows_count, columns_count;       


    SparseMatrix() : rows_count(0), columns_count(0) {}
    ~SparseMatrix() = default;
    SparseMatrix(size_t rows, size_t columns)
        : values(),
        row_indices(rows + 1),
        column_indices(),
        rows_count(rows),
        columns_count(columns) {}


    double& operator()(size_t row, size_t column)
    {
        size_t n1 = row_indices[row + 1];

        for (size_t k = row_indices[row]; k < n1; k++)
            if (column_indices[k] == column) return values[k];

        // adding a new element to a vector
        size_t k = row_indices[row + 1];
        values.insert(values.begin() + k, 0);
        column_indices.insert(column_indices.begin() + k, column);

        for (size_t i = row + 1, n2 = row_indices.size(); i < n2; i++)
            row_indices[i]++;
        return values[k];
    }

    double operator()(size_t row, size_t column) const
    {
        size_t n = row_indices[row + 1];
        for (size_t k = row_indices[row]; k < n; ++k)
            if (column_indices[k] == column)
                return values[k];
        return 0;
    }
};

template<class T>
int ILU_p(int n, int p, SparseMatrix<T> &U, SparseMatrix<T> &L)
{
    omp_set_num_threads(num_threads);
    std::vector <int> lev_u(n*n);
    std::vector <int> lev_l(n*n);
    std::vector <double> u_matr(n*n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j<n; j++) {
            u_matr[i*n + j] = U(i, j);
        }
    }

    #pragma omp parallel for collapse(2)
    for(int i =0; i < n; i++) {
        for (int j = 0;j < n; j++) {
            if (u_matr[i*n+j] != 0 or (i ==j)) {
                lev_u[i*n + j] = 0;
            } else{
                lev_u[i*n + j] = p+1;
            }
        }
    }


    for (int i = 1; i < n; i++)
    {
        for (int k = 0; k < i-1; k++){
            if (lev_u[i*n+k]<=p) {
                L(i, k) = U(i, k)/U(k,k);
                lev_l[i*n + k] = lev_u[i*n + k];
                for (int j = k; j < n; j++) {
                    U(i, j) -= L(i, k)*U(k, j);
                    lev_u[i*n + j] = lev_l[i*n + k] + lev_u[k*n + j] + 1;
                    if (lev_u[i*n + j] > p) {
                        U(i, j) = 0;
                    }
                }
            }
        }


    }

    return 0;
}


template<class T>
void printCSR(std::vector<size_t> row_ptr, std::vector<size_t> col_ind, std::vector<T> values, int n_rows) {
    int k = 0;
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_rows; j++) {
            if (k < col_ind.size() && col_ind[k] == j && k < row_ptr[i + 1]) {
                std::cout << values[k] << " ";
                k++;
            }
            else {
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;
    }
}


template<typename T>
int read_array(std::vector<T>& v, std::string str)
{
    int size;

    std::ifstream file(str);

    if (!file.is_open())
    {
        std::cout << "Error opening file" << std::endl;
        return 1;
    }

    file >> size;
    v.resize(size);

    for (int i = 0; i < size; i++){ file >> v[i]; }

    file.close();

    return 0;
}


int main()
{
    int n = 20;
    int p = 1;

    std::vector<size_t> rowPtr;
    std::vector<double> values;
    std::vector<size_t> colInd;

    read_array<size_t>(rowPtr, "indptr20.txt");
    read_array<double>(values, "data20.txt");
    read_array<size_t>(colInd, "indices20.txt");

    SparseMatrix<double> mat(n, n);
    //SparseMatrix<double> U(n, n);
    SparseMatrix<double> L(n, n);


    mat.values = values;
    mat.row_indices = rowPtr;
    mat.column_indices = colInd;
    mat.rows_count = n;
    mat.columns_count = n;

    auto start = std::chrono::steady_clock::now();

    ILU_p(n, p, mat, L);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;


    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    std::cout << "-------" << std::endl;
    printCSR(L.row_indices, L.column_indices, L.values, n);
    std::cout << "-------" << std::endl;
    printCSR(mat.row_indices, mat.column_indices, mat.values, n);
    std::cout << "-------" << std::endl;


    return 0;
}
