#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>

using namespace std;

extern "C" void dgehrd_(int* n, int ilo, int ihi, double* a, int* lda, double* tau, double* work, int* lwork, int* info);

int main(){


    int size = 1000;


    double* A = new double[size * size];

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            A[i * size + j] = i + j;
        }
    }

    int lwork = -1, info = 0, lda = size;
    double* tau = new double[size];
    double* exmp_work = new double[size];

    dgehrd_(&size, 1, size, A, &lda, tau, exmp_work, &lwork, &info);
    lwork = exmp_work[0];

    double* work = new double[lwork];


    auto start = chrono::steady_clock::now();

    sgehrd_(&size, 1, size, A, &lda, tau, work, &lwork, &info);

    auto end = chrono::steady_clock::now();
    cout << "Time in milliseconds:" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << "\n";
    cout << "Time in seconds:" << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << "\n";

    delete[] A;
    delete[] tau;
    delete[] work;
    delete[] exmp_work;

    return 0;
}
