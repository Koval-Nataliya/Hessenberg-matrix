#include <cmath>
#include <iostream>
#include <vector>
#include <time.h>
double sign(double x) {
    return x >= 0 ? 1 : -1;
}

double ** prod(double ** mass, double ** mass2, int n)
{
    auto **mass3 = new double *[n];
    for (int i = 0; i < n; i++) {
        mass3[i] = new double[n];
        for (int j = 0; j < n; j++) {
            mass3[i][j] = 0;
            for (int k = 0; k < n; k++) {
                mass3[i][j] += mass[i][k] * mass2[k][j];
            }
        }
    }
    return mass3;
}

int main () {
    int i, j, m;
    double s, r;
    m = 512;
    auto ** a = new double*[m];

    for(i = 0; i < m; i++) {
        a[i] = new double[m];
        for (j = 0; j < m; j++) {
            a[i][j] = rand();
        }
    }
    unsigned int start_time = clock();

    for (i = 0; i < m - 2; i++) {
        double norm = 0;
        for (j = i + 1; j < m; j++) {
            norm += a[j][i] * a[j][i];
        }
        norm = sqrt(norm);
        s = sign(a[i + 1][i]) * norm;
        r = sqrt(2 * a[i + 1][i] * s + 2 * s * s);

        std::vector<double > w(m, 0);

        w[i + 1] = (a[i + 1][i] + s) * (1 / r);
        for (j = i + 2; j < m; j++) {
            w[j] = a[j][i] * (1 / r);
        }

        auto ** eye = new double*[m]; //единичная матрица
        for(j = 0; j < m; j++) {
            eye[j] = new double[m];
            for (int k = 0; k < m; k++) {
                if (j == k)
                    eye[j][k] = 1;
                else
                    eye[j][k] = 0;
            }
        }

        for(int k = 0; k < m; k++) {
            for(int t = 0; t < m; t++) {
                eye[k][t] -= 2 * w[k] * w[t]; // матрица h
            }
        }

        a = prod(eye, a, m);
        a = prod(a, eye, m);
    }

    unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time;
    /* print the matrix A after calling hessenberg */
//    printf("A = \n");
//    for(i = 0; i < m; i++) {
//        for(j = 0; j < m; j++) {
//            double sign = std::copysignf(1, a[i][j]);
//            a[i][j] = sign * round(std::copysignf(1, a[i][j]) * a[i][j] * 100) / 100;
//            std::cout<<a[i][j]<<' ';
//        }
//        printf("\n");
//    }
//    printf("\n");
    std::cout << search_time / CLOCKS_PER_SEC << std::endl;

//    /* free memory */
    for(i = 0; i < m - 1; i++) {
        delete[] a[i];
    }
    delete[] a;
    return 0;
}
