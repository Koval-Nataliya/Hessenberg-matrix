#include<iostream>
#include<vector>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include <ctime>
#include <chrono>
#include <iomanip>


// создание пустой матрицы
std::vector<float> create_matrix(int size)
{
    std::vector<float> a(size * size);
    for (int i = 0; i < size * size; i ++) {
        a[i] = rand();
    }
    return a;
}

int dot(std::vector<float>& a, std::vector<float>& b, std::vector<float>& ans, int size, float k)
{

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            ans[i * size + j] = k * a[i] * b[j];
        }
    }
    return 0;
}

// умножение строки на матрицу
int dot_matrix_vector1(std::vector<float>& matrix, std::vector<float>& b, std::vector<float>& ans, int s)
{
    for (int k = 0; k < s; k++)
    {
        ans[k] = 0;
        for (int i = 0; i < s; i++) {
            ans[k] += matrix[i * s + k] * b[i];
        }
    }
    return 0;
}

// умножение матрциы на столюец
int dot_matrix_vector2(std::vector<float>& matrix, std::vector<float>& b, std::vector<float>& ans, int s)
{
    for (int k = 0; k < s; k++)
    {
        ans[k] = 0;
        for (int i = 0; i < s; i++) {
            ans[k] += matrix[i + k * s] * b[i];
        }
    }
    return 0;
}
void is_unitary(float ** arr, int size) {
    double eps = 1e-7;
    for (int i =0; i < size; i++) {
        for (int j =0; j < size; j++) {
            if ((arr[i][j] - arr[j][i]) > eps) {
                std::cout << "Matrix is not unitary!" << '\n';
                return;
            }
        }
    }
    std::cout << "Matrix is unitary!" << '\n';
}
// преобразование матрицы по алгоритму Хаусхолдера
int hausholders_transformation(std::vector<float> &matrix, int size)
{
    for (int i = 0; i < size - 2; i++)
    {
        float s = 0;
        for (int j = i + 2; j < size; j++) {
            s += fabs(matrix[j * size + i]) * fabs(matrix[j * size + i]);
        }
        float norm_1 = sqrtf(fabs(matrix[i + (i + 1) * size]) * fabs(matrix[i + (i + 1) * size]) + s);
        std::vector<float> x(size);
        for (int j = 0; j < size; j++) {
            if (j < i + 1) {
                x[j] = 0;
            } else if(j == i + 1){
                x[j] = matrix[i + j * size] - norm_1;
            } else {
                x[j] = matrix[i + j * size];
            }
        }

        float norm_x = sqrtf(fabs(x[i + 1]) * fabs(x[i + 1]) + s);
        for(int j = 0; j < size; j++) {
            x[j] = x[j] / norm_x;
        }

//        check unitar
//        Проверяем каждую матрицу U на унитарность, тк произведение унитраных матриц унитарно
//        std::vector<float> check(size * size);
//        dot(x, x, check, size, 2);
//
//        auto ** u = new float*[size];
//
//        for(int t = 0; t < size; t++) {
//            u[t] = new float[size];
//            for (int j = 0; j < size; j++) {
//                if (t == j) {
//                    u[t][j] = 1 - check[t * size + j];
//                } else {
//                    u[t][j] = 0 - check[t * size + j];
//                }
//            }
//        }
//        is_unitary(u, size);

        std::vector<float> y(size);
        dot_matrix_vector1(matrix, x, y, size);

        std::vector<float> z(size);
        dot_matrix_vector2(matrix, x, z, size);

        std::vector<float> tmp(size * size);
        dot(x, y, tmp, size, 2);

        std::vector<float> tmp1(size * size);
        dot(z, x, tmp1, size, 2);

        float b = 0;
        for(int j = 0; j < size; j++) {
            b += y[j] * x[j];
        }
        b *= 4;

        std::vector<float> prom(size * size);
        dot(x, x, prom, size, b);

        for(int j = 0; j < size * size; j++) {
            matrix[j] = matrix[j] - tmp1[j] - tmp[j] + prom[j];
        }
    }
    return 0;
}


void is_hessenberg(std::vector<float> &matrix, int size) {
    double eps = 1e-7;
    for(int i = 2; i < size; i++) {
        for(int j = 0; j < i - 1; j++) {
            if (round(fabs(matrix[i * size + j]) * 100)/100 > eps) {
                std::cout << "Matrix is not Hessenberg!" << '\n';
                return;
            }
        }
    }
    std::cout << "Matrix is Hessenberg!" << '\n';
}

int main()
{
    int matrix_size= 2048;

    std::vector<float> matrix(matrix_size * matrix_size);
    for (int i = 0; i < matrix_size * matrix_size; i++) {
        matrix[i] = rand();
    }
//    std::vector<float> matrix = {-602,  -84,  681,  867,  613, -298, -947,   70, -232,  616,  582, 337, -932, -584, -673,   83,  143,  564, -931,  909, -419,   84,
//                                     -236, -773, -419, -847,  745,  382, -967,  852,  878,  693, -490, 898, -998,   15, -700,  606, -147, -663, -900,  179, -955,  -82,
//                                     475,  765, -457, -420, -425,  120, -831,  246, -212,  781, -105, 347, -179, -904,  775, -523,   42,  732,  583, -507,  952,  443,
//                                     960,  150, -981, -272, -851, -351,  -74,  436,   -5, -660, -503, 200, -958,  261,  116,  -94, -732, -183,  416,  819,  936, -878,
//                                     204, -415,  134,  -68,  -10, -558, -875,  812,  344, -273, -785, 186, -681, -411,  337,  815,  194, -213, -128,  440, -554, -264,
//                                     944,  423, -627,  527,  934, -619,  375,  248,  683,  908, -398, -419, -545, -138, -806, -692, -239,  847, -270,  154, -517,  283,
//                                     -578, -107, -714,    6, -519, -252,  270, -677,  547,  437,  276, -467, -428, -102,  647, -553, -141,  892,  368,  847, -284, -755,
//                                     995,  829,  641, -917,  325, -853, -702,  178,  272,  105, -147, -658, -725,  427,  883,  567,  403, -389,  673,  617,  590, -805,
//                                     -91,  693,  253,  -30, -974,  407,  314, -805, -178,  964,  662, -460, -639,    0, -385, -629,  190, -493, -982, -114, -710, -808,
//                                     -429,  678, -439,  268,   45, -450,    8,  -50,  -96,  602, -458, 574, -628,  178, -117, -572,  761, -870,  304, -592, -421,   60,
//                                     843,  505, -785,  413,  990, -395,  673, -365,  114, -779, -658, 18, -621,  361,  710,  932, -787,  670, -541, -736,  525,  153,
//                                     -184,  906,  425,  518,  -94, -822,  456, -405, -466, -232,  811, -965,   69,   44};
   // std::vector<float> matrix = {4, 2, 23, 3, 1, -4, 12, 3, 8, 2, 2, 9, 1, -4, 7, -3, 5, -5, 11, -8, 1, 3, 1, 2, 1, 3, 1, 1, 2, 3, 21, 43, 5, -3, 6, 11};

    auto start = std::chrono::steady_clock::now();

    hausholders_transformation(matrix, matrix_size);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;


    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
//    for (int i = 0; i < matrix_size * matrix_size; i++) {
//        if (i % matrix_size == 0) {
//            std::cout<<'\n';
//        }
//        std::cout << roundf(matrix[i] * 100) / 100 << ' ';
//    }
//    std::cout<<'\n';
    //is_hessenberg(matrix, matrix_size);

    return 0;
}

