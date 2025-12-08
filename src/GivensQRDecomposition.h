#pragma once

#include <iostream>
#include <cmath>

#include "TMatrix.h"

template <typename T>
void Givens_QR_decomposition(const TMatrix<T> &A, TMatrix<T> &Q, TMatrix<T> &R)
{
    auto n = A.Size(), m = A[0].Size();
    R = A;
    bool isFirst = true;
    for (auto j = 0; j < m - 1; j++)
    { // зануляем весь столбец j
        for (auto i = j + 1; i < n; i++)
        { // c j+1 строк
            std::cout << "i=" << i << ", j=" << j << std::endl;

            auto sqrt = std::sqrt(R[j][j] * R[j][j] + R[i][j] * R[i][j]);
            double c = R[j][j] / sqrt, s = -R[i][j] / sqrt;
            std::cout << "c=" << c << ", s=" << s << std::endl;

            for (auto k = j; k < m; k++)
            { // вращаем только две строки -- i и j
                auto temp = R[j][k] * c - R[i][k] * s;
                R[i][k] = R[j][k] * s + R[i][k] * c;
                R[j][k] = temp;
            }
            std::cout << "R:" << R;
            std::cout << std::endl;

            if (isFirst)
            {
                Q = TMatrix<T>(n, n, 0);
                for (auto i = 0; i < n; i++)
                {
                    Q[i][i] = 1;
                }
                Q[0][0] = Q[1][1] = c;
                Q[0][1] = -s;
                Q[1][0] = s;
                isFirst = false;
                // std::cout << "Q:" << Q;
                // std::cout << std::endl;
            }
            else
            {
                for (auto k = j; k < n; k++)
                { // меняются только две строки -- i и j
                    auto temp = Q[j][k] * c - Q[i][k] * s;
                    Q[i][k] = Q[j][k] * s + Q[i][k] * c;
                    Q[j][k] = temp;
                }
            }
            std::cout << "Q:" << Q;
            std::cout << std::endl;
            std::cout << "-------------------\n";
        }
    }

    // транспонировать Q
    for (auto i = 0; i < n - 1; i++)
    {
        for (auto j = i + 1; j < m; j++)
        {
            std::swap(Q[i][j], Q[j][i]);
        }
    }
    // std::cout << "Q:" << Q;
    // std::cout << std::endl;
}
