#include <iostream>
#include <vector>
#include <algorithm>    // std::find
#include <omp.h>
#include <chrono>
using namespace std;
int main() {
    //три вектора со значениями и индексами
    vector<int>  rows, cols;
    vector<double> values;

    vector<int>  rows_L, cols_L;
    vector<double> values_L;

    vector<int>  rows_U, cols_U;
    vector<double> values_U;

    double sum = 0;
    //fill values
    int n = 100;
    int num_fill = 40;
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

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (i == j) {
                rows_L.push_back(i);
                cols_L.push_back(i);
                values_L.push_back(1);

                rows_U.push_back(i);
                cols_U.push_back(i);
                values_U.push_back(1);
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
            if (it != it2+1) {
                for(int k = 0; k <j-1; k++) {
                    auto iter_tmp_L1 = find(rows_L.begin(), rows_L.end(), i);
                    auto iter_tmp_L2 = find(cols_L.begin(), cols_L.end(), k);

                    auto iter_tmp_U1 = find(rows_U.begin(), rows_U.end(), k);
                    auto iter_tmp_U2 = find(cols_U.begin(), cols_U.end(), j);
                    int distL1 = distance(rows_L.begin(), iter_tmp_L1);
                    int distL2 = distance(cols_L.begin(), iter_tmp_L2);
                    int distU1 = distance(rows_U.begin(), iter_tmp_U1);
                    int distU2 = distance(cols_U.begin(), iter_tmp_U2);
                    if (distL1 == distL2 && distU1 == distU2) {
                        sum += values_L[distL1] * values_U[distU1];
                    }
                }
                int ind = distance(cols.begin(), it);
                auto iter_U = find(rows_U.begin(), rows_U.end(), j);
                int d = distance(rows_U.begin(), iter_U);
                values_L.push_back((values[ind] - sum)/values_U[d]);
                rows_L.push_back(i);
                cols_L.push_back(j);
            }
            sum = 0;
        }
        for(int j = i+1; j < n; j++) {
            auto it = find(it1, it2, j);
            if (it != it2+1) {
                for(int k = 0; k <j-1; k++) {
                    auto iter_tmp_L1 = find(rows_L.begin(), rows_L.end(), i);
                    auto iter_tmp_L2 = find(cols_L.begin(), cols_L.end(), k);

                    auto iter_tmp_U1 = find(rows_U.begin(), rows_U.end(), k);
                    auto iter_tmp_U2 = find(cols_U.begin(), cols_U.end(), j);
                    int distL1 = distance(rows_L.begin(), iter_tmp_L1);
                    int distL2 = distance(cols_L.begin(), iter_tmp_L2);
                    int distU1 = distance(rows_U.begin(), iter_tmp_U1);
                    int distU2 = distance(cols_U.begin(), iter_tmp_U2);
                    if (distL1 == distL2 && distU1 == distU2) {
                        sum += values_L[distL1] * values_U[distU1];
                    }
                }
                int ind = distance(cols.begin(), it);
                values_U.push_back(values[ind] - sum);
                rows_U.push_back(i);
                cols_U.push_back(j);
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
