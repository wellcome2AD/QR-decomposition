#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "QRDecomposition.h"

int main1() {
	const int N = 3;
	TMatrix<float> A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} };
	const TMatrix<float> A_copy = A;
	const TVector<float> b = { 1, 8, 5 };
	TVector<float> x = QR_decomposition(A, b);
	std::cout << "x\n" << x;

	// проверка
	const float eps = 0.000001;
	TVector<float> A_x = A_copy * x;
	for (auto i = 0; i < A_x.Size(); ++i) {
		float bi = b[i];
		if (abs(bi - A_x[i]) >= eps) {
			printf("error in %I32d: %f != %f\n", i, A_x[i], b[i]);
		}
	}

	return 0;
}
