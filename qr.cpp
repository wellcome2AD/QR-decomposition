#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "invertible_matrix.h"

typedef std::vector<float> fvector;
typedef std::vector<std::vector<float>> fmatrix;

float sign(float x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 1;
}

fvector vec_mult(fvector v, float a) {
    fvector res(v.size());
    for(ptrdiff_t i = 0; i < v.size(); ++i) {
        res[i] = v[i] * a;
    }
    return res;
}

fmatrix vec_mult(fvector v1, fvector v2) {
    ptrdiff_t N = v1.size();
    fmatrix res(N, fvector(N));
    for(ptrdiff_t i = 0; i < N; ++i) {
        for(ptrdiff_t j = 0; j < N; ++j) {
            res[i][j] = v1[i] * v2[j];
        }
    }
    return res;
}

fmatrix matr_diff(fmatrix m1, fmatrix m2) {
    ptrdiff_t N = m1.size();
    for(ptrdiff_t i = 0; i < N; ++i) {
        for(ptrdiff_t j = 0; j < N; ++j) {
            m1[i][j] -= m2[i][j];
        }
    }
    return m1;
}

fmatrix matr_mult(fmatrix m1, fmatrix m2) {
    ptrdiff_t N = m1.size();
    ptrdiff_t M = m2[0].size();
    fmatrix res(N, fvector(M, 0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            for (int k = 0; k < N; ++k) {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return res;
}

fmatrix count_H(float beta, float mu, fvector w) {
    ptrdiff_t N = w.size();
    fmatrix E(N, fvector(N, 0));
    for(ptrdiff_t i = 0; i < N; ++i) {
        E[i][i] = 1;
    }
    return matr_diff(E, vec_mult(vec_mult(w, 2), w));
}

void print_matr(fmatrix m) {
    ptrdiff_t N = m.size();
    ptrdiff_t M = m[0].size();
    for(ptrdiff_t i = 0; i < N; ++i) {
        std::cout << std::left << std::setw(15);
        for(ptrdiff_t j = 0; j < M; ++j) {
            std::cout << m[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main()
{
    const ptrdiff_t N = 3;
    fmatrix A = {{1, -2, 1}, {2, 0, -3}, {2, -1, -1}};
    const fmatrix A_copy = A;
    const fmatrix b = { {1}, {8}, {5} };
    fmatrix x = { {0}, {0}, {0} };
    fmatrix Q;
    bool isFirst = true;
    for(ptrdiff_t k = 0; k < N - 1; ++k) { // k-ый шаг алгоритма
        fvector A_k_col(N);
        float sum_by_k_col = 0;
        for(ptrdiff_t str = k; str < N; ++str) { // нужен не весь столбец, а только начиная с k-ой строки
            A_k_col[str] = A[str][k];
            sum_by_k_col += A[str][k]*A[str][k];
        }

        float beta = sign(-A[k][k]) * sqrt(sum_by_k_col);
        printf("beta%ld = %e\n", k, beta);
        float mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * A[k][k]);
        printf("mu%ld = %e\n", k, mu);

        A_k_col[k] -= beta;
        fvector w = vec_mult(A_k_col, mu);
        for(ptrdiff_t index = 0; index < k; ++index) { // первые k позиций - нулевые
            w[index] = 0;
        }

        fmatrix H = count_H(beta, mu, w);
        printf("H%ld\n", k+1);
        print_matr(H);
        std::cout << std::endl;
        
        if (isFirst) {
            isFirst = false;
            Q = H;
        } else {
            Q = matr_mult(Q, H);
        }

        // умножить Hk на A, получим новую Ak
        A = matr_mult(H, A);
        printf("A%ld\n", k+1);
        print_matr(A);
        std::cout << std::endl;
    }
    
    printf("Q\n");
    print_matr(Q);
    std::cout << std::endl;

    // обратный ход
    fmatrix Q_invert = Mreverse(Q);
    printf("Q_invert\n");
    print_matr(Q_invert);
    std::cout << std::endl;

    fmatrix Q_invert_b = matr_mult(Q_invert, b);
    printf("Q_invert_b\n");
    print_matr(Q_invert_b);
    std::cout << std::endl;

    x[N - 1][0] = Q_invert_b[N - 1][0] / A[N - 1][N - 1];
    for (ptrdiff_t i = N - 2; i >= 0; --i) {
        float sum = 0;
        for (ptrdiff_t j = i + 1; j < N; ++j) {
            sum += A[i][j] * x[j][0];
        }
        x[i][0] = (Q_invert_b[i][0] - sum) / A[i][i];
    }

    printf("x\n");
    print_matr(x);
    std::cout << std::endl;

    // проверка
    const float eps = 0.000001;
    fmatrix A_x = matr_mult(A_copy, x);
    for (ptrdiff_t i = 0; i < A_x.size(); ++i) {
        for (ptrdiff_t j = 0; j < A_x[0].size(); ++j) {
            if (abs(b[i][j] - A_x[i][j]) >= eps) {
                printf("error in %ld %ld: %f != %f\n", i, j, A_x[i][j], b[i][j]);
            }
        }
    }

    return 0;
}
