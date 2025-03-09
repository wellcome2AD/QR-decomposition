#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

typedef std::vector<float> fvector;
typedef std::vector<std::vector<float>> fmatrix;

float sign(float x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 1;
}

fvector vec_mult(fvector v, float a) {
    fvector res(v.size());
    for(size_t i = 0; i < v.size(); ++i) {
        res[i] = v[i] * a;
    }
    return res;
}

fmatrix vec_mult(fvector v1, fvector v2) {
    size_t N = v1.size();
    fmatrix res(N, fvector(N));
    for(size_t i = 0; i < N; ++i) {
        for(size_t j = 0; j < N; ++j) {
            res[i][j] = v1[i] * v2[j];
        }
    }
    return res;
}

fmatrix matr_diff(fmatrix m1, fmatrix m2) {
    size_t N = m1.size();
    for(size_t i = 0; i < N; ++i) {
        for(size_t j = 0; j < N; ++j) {
            m1[i][j] -= m2[i][j];
        }
    }
    return m1;
}

fmatrix matr_mult(fmatrix m1, fmatrix m2) {
    size_t N = m1.size();
    fmatrix res(N, fvector(N, 0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return res;
}

fmatrix count_H(float beta, float mu, fvector w) {
    size_t N = w.size();
    fmatrix E(N, fvector(N, 0));
    for(size_t i = 0; i < N; ++i) {
        E[i][i] = 1;
    }
    return matr_diff(E, vec_mult(vec_mult(w, 2), w));
}

void print_matr(fmatrix m) {
    size_t N = m.size();
    for(size_t i = 0; i < N; ++i) {
        std::cout << std::left << std::setw(15);
        for(size_t j = 0; j < N; ++j) {
            std::cout << m[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main()
{
    const size_t N = 3;
    fmatrix A = {{1, -2, 1}, {2, 0, -3}, {2, -1, -1}};
    const fmatrix A_copy = A;
    const fvector b = {1, 8, 5}, x = {0, 0, 0};
    fmatrix Q;
    bool isFirst = true;
    for(size_t k = 0; k < N - 1; ++k) { // k-ый шаг алгоритма
        fvector A_k_col(N);
        float sum_by_k_col = 0;
        for(size_t str = k; str < N; ++str) { // нужен не весь столбец, а только начиная с k-ой строки
            A_k_col[str] = A[str][k];
            sum_by_k_col += A[str][k]*A[str][k];
        }

        float beta = sign(-A[k][k]) * sqrt(sum_by_k_col);
        printf("beta%ld = %e\n", k, beta);
        float mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * A[k][k]);
        printf("mu%ld = %e\n", k, mu);

        A_k_col[k] -= beta;
        fvector w = vec_mult(A_k_col, mu);
        for(size_t index = 0; index < k; ++index) { // первые k позиций - нулевые
            w[index] = 0;
        }

        fmatrix H = count_H(beta, mu, w);
        printf("H[%ld]\n", k);
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
        printf("A[%ld]\n", k);
        print_matr(A);
        std::cout << std::endl;
    }
    
    printf("Q\n");
    print_matr(Q);
    std::cout << std::endl;
    
    return 0;
}
