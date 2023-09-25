#include <iostream>
#include <vector>
#include <algorithm>    // std::find
#include <omp.h>
#include <chrono>
using namespace std;
int num_threads = 4;
int main() {
    //три вектора со значениями и индексами
    omp_set_num_threads(num_threads);
    vector<int>  rows, cols;
    vector<double> values;
    double sum = 0;
    //fill values
    int n = 500;
    int num_fill = 60;
    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < num_fill; i++) {
        double v = rand() % 100;
        int r = rand() % n;
        int c = rand() % n;
        values.push_back(v);
        rows.push_back(r);
        cols.push_back(c);
    }

    sort(rows.begin(), rows.end());
    sort(cols.begin(), cols.end());
//    for(int i = 0; i<10;i++){
//        cout<<rows[i]<<' ';
//    }
    double L[n][n], U[n][n];
    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
            if (i == j) {
                L[i][j] = 1;
                U[i][j] = 1;
            }
        }
    }
    auto it1 = cols.begin();
    auto it2 = cols.begin();
    auto start = std::chrono::steady_clock::now();

    for(int i = 0; i < n; i++) {
        int cnt = count(rows.begin(), rows.end(), i);
        it2 = it2 + cnt;
        sum = 0;
        for(int j = 0; j<=i-1; j++) {
            auto it = find(it1, it2, j);
            if (it == it2+1) {
                L[i][j] = 0;
            } else {
                #pragma omp parallel for num_threads(num_threads)
                for(int k = 0; k <j-1; k++) {
                    sum += L[i][k]*U[k][j];
                }
                int ind = distance(cols.begin(), it);
                L[i][j] = (values[ind] - sum)/U[j][j];
            }
            sum = 0;
        }
        L[i][i] = 1;
        for(int j = i+1; j < n; j++) {
            auto it = find(it1, it2, j);
            if (it == it2+1) {
                U[i][j] = 0;
            } else {
                #pragma omp parallel for num_threads(num_threads)
                for(int k = 0; k <j-1; k++) {
                    sum += L[i][k]*U[k][j];
                }
                int ind = distance(cols.begin(), it);
                U[i][j] = (values[ind] - sum);
            }
            sum = 0;
        }
        it1 = it2;
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;


    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    
    return 0;
}
