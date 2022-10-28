#include <cmath>
#include <iostream>
#include <vector>
#include <time.h>
#include <cstring>

const int b = 64;
const int n = 1024;

void multiply_block(double *X, double *Y, double *res) {
    for (int j = 0; j < b; ++j)
        for (int s = 0; s < b; ++s)
            for (int i = 0; i < b; ++i)
                res[i + j * b] += X[i + s * b] * Y[s + j * n];
}

void copy_a(double *block_a, double *A) {
    for (int j = 0; j < b; ++j)
        for (int i = 0; i < b; ++i)
            block_a[i + j * b] = A[i + j * n];
}

void copy_c(double *block_c, double *C) {
    for (int j = 0; j < b; ++j)
        for (int i = 0; i < b; ++i)
            C[i + j * n] = block_c[i + j * b];
}

void matmul(double *A, double *B, double *C) {

    static double block_c[b * b];
    static double block_a[b * b];

    for (int j = 0; j < n; j += b) {
        for (int i = 0; i < n; i += b) {
            memset(block_c, 0, b * b * sizeof(double));
            for (int k = 0; k < n; k += b) {
                copy_a(block_a, &A[i + k * n]);
                multiply_block(block_a, &B[k + j * n], block_c);
            }
            copy_c(block_c, &C[i + j * n]);
        }
    }
}
//double* transpose(double* h){
//    int t;
//    double *a = (double *)malloc(n * n * sizeof(double));
//
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < n; ++j) {
//            a[i*n + j] = h[i +j *n];
//        }
//    }
//    return a;
//}

int main () {
    int i, j, m;
    double s, r;
    m = n;
    double *a = (double *)malloc(n * n * sizeof(double));
    double *ans = (double *)malloc(n * n * sizeof(double));
    for(i = 0; i < m * m; i++) {
        a[i] = rand();
    }
    unsigned int start_time = clock();

    for (i = 0; i < m - 2; i++) {
        double norm = 0;
        for (j = i + 1; j < m; j++) {
            norm += a[m * j + i] * a[m * j + i];
        }
        norm = sqrt(norm);
        s = std::copysignf(1, a[(i + 1)*m + i]) * norm;
        r = sqrt(2 * a[(i + 1)*m + i] * s + 2 * s * s);

        std::vector<double > w(m, 0);

        w[i + 1] = (a[(i + 1) * m + i] + s) * (1 / r);
        for (j = i + 2; j < m; j++) {
            w[j] = a[j*m + i] * (1 / r);
        }

        double *eye= (double *)malloc(n * n * sizeof(double));

        for(j = 0; j < m * m; j++) {
            if (j % (m + 1) == 0) {
                eye[j] = 1;
            }
            else {
                eye[j] = 0;
            }
        }
        for(int k = 0; k < m; k++) {
            for(int t = 0; t < m; t++) {
                eye[k * m + t] -= 2 * w[k] * w[t]; // матрица h
            }
        }
        matmul(a, eye, ans);
        a = ans;
        matmul(eye, a, ans);
        a = ans;
    }

    unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time;
    std::cout << search_time / CLOCKS_PER_SEC << std::endl;
    return 0;
}
