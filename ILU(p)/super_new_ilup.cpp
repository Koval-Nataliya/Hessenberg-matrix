#include <iostream>
#include <vector>
#include <algorithm>    // std::find
#include <utility>
#include <map>
using namespace std;

double find_pair(vector<pair<double, pair<int, int>>> a, int k) {
    for(int i = 0; i < a[0].first; i++) {
        if (k == a[i].second.first) {
            return a[i].first;
        }
    }
    return 0;
}

int main() {
    int p = 2;
    int n = 6;
    //массив разреженных строк, где первое значение - это значение, и пара место-уровень
    std::vector<pair<double, pair<int, int>>> arr[n];
    std::vector<pair<double, pair<int, int>>> U[n];
    std::vector<pair<double, int>> L[n];

    //length of vector;
    int r = 0;
    //первый элемент - 0, (количество, -1)
    for (int i = 0; i < n; i++) {
        r = rand() % 5;
        arr[i].push_back(make_pair(0 ,make_pair(r, -1)));
    }
    for (int i = 0; i<n; i++) {
        for(int j = 0; j <= arr[i][0].first; j++) {
            double num = rand() % 100;
            r = rand() % 6;
            arr[i].push_back(make_pair(num, make_pair(r, 1)));
        }
    }
    //fill levels
    for(int i =0; i < n;i++) {
        for(int j = 0; j<arr[i][0].second.second; j++) {
            if (i != j)
                arr[i][j].second.second = p+1;
        }
    }
    //fill U as arr
    for (int i = 0; i < n; i++) {
        for(int j = 0; j < arr[i][0].first; j++) {
            if (j == 0) {
                U[i][0].first = 0;
                U[i][0].second.first = arr[i][0].second.first;
                U[i][0].second.second = -1;
            } else {
                U[i][j].first = arr[i][j].first;
                U[i][j].second.first = arr[i][j].second.first;
                U[i][j].second.second = arr[i][j].second.second;
            }
        }
    }
    for (int i =0; i < n- 1; i++) {
        for (int j = i+1; j < n; j++) {
            double tmp = find_pair(arr[j], i);
            double tmp1 = find_pair(U[j], i);
            if (tmp1 != 0 && tmp !=0) {
                double koef = tmp / tmp1;
                for (int k = i; k < n; k++){
                    double t = find_pair(U[j], k);
                    if (t != 0) {
                        U[j][k].first = 0;
                        U[j][k].second.first = 0;
                        U[j][k].second.second = 0;
                    }
                    L[i].push_back(make_pair(koef, j));
                }
            }

        }
    }
    return 0;
}
