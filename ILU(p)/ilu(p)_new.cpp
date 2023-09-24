#include <iostream>
#include <vector>
#include <algorithm>    // std::find
using namespace std;
int main() {
    //три вектора со значениями и индексами
    vector<int>  rows, cols;
    vector<double> values;
    double sum = 0;
    //fill values
    int n = 6;
    for(int i = 0; i < 10; i++) {
        double v = rand() % 100;
        int r = rand() % 6;
        int c = rand() % 6;
        values.push_back(v);
        rows.push_back(r);
        cols.push_back(c);
        cout<<values[i]<<' ';
    }

    sort(rows.begin(), rows.end());
    sort(cols.begin(), cols.end());
    for(int i = 0; i<10;i++){
        cout<<rows[i]<<' ';
    }
    double L[n][n], U[n][n];
    for(int i; i < n; i++) {
        for(int j = 0; j < n; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
            if (i == j) {
                L[i][j] = 1;
            }
        }
    }
    auto it1 = cols.begin();
    auto it2 = cols.begin();
    for(int i = 0; i < n; i++) {
        int cnt = count(rows.begin(), rows.end(), i);
        it2 = it2 + cnt;
        sum = 0;
        for(int j = 0; j<=i-1; j++) {
            auto it = find(it1, it2, j);
            if (it == it2+1) {
                L[i][j] = 0;
            } else {
                for(int k = 0; k <j-1; k++) {
                    sum += L[i][k]*U[k][j];
                }
                int ind = distance(cols.begin(), it);
                L[i][j] = (values[ind] - sum)/U[j][j];
            }
        }
        L[i][i] = 1;
        for(int j = i+1; j < n; j++) {
            auto it = find(it1, it2, j);
            if (it == it2+1) {
                U[i][j] = 0;
            } else {
                for(int k = 0; k <j-1; k++) {
                    sum += L[i][k]*U[k][j];
                }
                int ind = distance(cols.begin(), it);
                U[i][j] = (values[ind] - sum);
            }
        }
        it1 = it2;
    }

    for(int i = 0; i<n; i++) {
        cout<<'\n';
        for(int j =0; j< n; j++) {
            cout<<
            L[i][j]<<' ';
        }
    }
    return 0;
}
