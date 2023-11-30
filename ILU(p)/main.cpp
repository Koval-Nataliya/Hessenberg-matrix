#include <iostream>
#include <vector>
#include <algorithm>    // std::find
#include <utility>
#include <map>
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
    int number_of_str = 10;
    struct SparseVec vec[number_of_str];
    int n = 0;
    //filling matrix;
    for(int i = 0; i<number_of_str; i++) {
        n = rand()%5;
        vec[i].nz = n;
        //чтобы был диагональный элемент
        vec[i].ind.push_back(i);
        vec[i].val.push_back(rand()%100);
        vec[i].levels.push_back(1);
        //остальные элементы
        for (int j = 0; j < n-1; j++) {
            vec[i].ind.push_back(vec[i].nz - int((j - i)/2) + rand()%4);
            vec[i].val.push_back(rand()%100);
            vec[i].levels.push_back(1);
        }
        std::sort(vec[i].ind.begin(), vec[i].ind.end());
    }

    //main cycle
    for(int i = 0; i< number_of_str; i++) {
        for (int j = 0; j <number_of_str;j++) {
            vec[j] = axy(vec[i], vec[j], i);
        }
    }

//    for (int i = 0; i < number_of_str; i++) {
//        std::cout<<"nz "<<vec[i].nz << " ind ";
//        for (int j = 0; j < vec[i].nz; j++) {
//            std::cout<<vec[i].ind[j]<<' ';
//        }
//        std::cout<<"val: ";
//        for (int j = 0; j < vec[i].nz; j++) {
//            std::cout<<vec[i].val[j]<<' ';
//        }
//        std::cout<<"lev ";
//        for (int j = 0; j < vec[i].nz; j++) {
//            std::cout<<vec[i].levels[j]<<' ';
//        }
//        std::cout<<'\n';
//    }
    return 0;
}
