#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <map>
#include <vector>
#include "IQRSolver.h"

#include "Householder_method/HouseholderBasic.h"
#include "Householder_method/HouseholderWithoutMatrixMul.h"
#include "Householder_method/HouseholderWithNormW.h"

#include "Givens_method/GivensBasic.h"
#include "Givens_Method/GivensQRInOneMatrix.h"
#include "Givens_Method/GivensSIMD.h"

#include "matrix.h"

typedef double currentType;
int thread_num = 8;

template <typename T>
void printResultWithExpected(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& Q, const std::vector<std::vector<T>>& R)
{
	// проверка на соответствие дл€ малых тестов с помощью базовой реализации
	IQRSolver<currentType>* basicMethod = new HouseholderMethodBasic<currentType>();
	std::vector<std::vector<T>> expecQ, expecR;
	basicMethod->QR_decomposition(A, expecQ, expecR);
	std::cout << "Q:";
	printMatrix(Q);
	std::cout << std::endl;
	std::cout << "expected Q:";
	printMatrix(expecQ);
	std::cout << std::endl;

	std::cout << "R:";
	printMatrix(R);
	std::cout << std::endl;
	std::cout << "expected R:";
	printMatrix(expecR);
	std::cout << std::endl;

	std::cout << "Q * R:";
	printMatrix(multiplyMatrix(Q, R));
	std::cout << std::endl;
	std::cout << "A:";
	printMatrix(A);
	std::cout << std::endl;
}

template <typename T>
void QR_decomposition_test_with_generated_matrix(IQRSolver<T>* solver, int N, bool writeToFile)
{
	std::cout << "test size: " << N << std::endl;
	std::vector<std::vector<T>> A = generateMatrix<T>(N, N);
	std::vector<std::vector<T>> Q(N, std::vector<T>(N, 0)), R(N, std::vector<T>(N, 0));
	// std::vector < std::vector <currentType>> A = { {1,2,3}, {4, 5, 6}, {7, 8, 9} }, R, Q;

	double start = omp_get_wtime();
	solver->QR_decomposition(A, Q, R);
	double end = omp_get_wtime();

	std::cout << "time: " << end - start << std::endl;
	std::cout << "abs error: " << Fnorm<T>(substractMatrix<T>(multiplyMatrix<T>(Q, R), A)) << std::endl;
	std::cout << "rel error: " << Fnorm<T>(substractMatrix<T>(multiplyMatrix<T>(Q, R), A)) / Fnorm<T>(A) << std::endl << std::endl;

	// printResultWithExpected(A, Q, R);

	if (writeToFile) {
		writeMatrixToFile<T>("test_data\\matrixA_" + std::to_string(N) + ".txt", A);
		writeMatrixToFile<T>("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile<T>("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
}

struct testParams
{
	IQRSolver<currentType>* solver;
	std::string description;
	std::vector<int> sizes;
};

void QRtests()
{
	auto methods = std::map<int, testParams>{};

	/*methods[0] = {
		new HouseholderMethodBasic<currentType>(),
		"Householder basic version",
		{ 100, 200, 300, 400, 500 },
	};
	
	methods[1] = {
		new HouseholderMethodWithoutMatrixMults<currentType>(),
		"Householder without matrix multiplications",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};
	
	methods[2] = {
		new HouseholderMethodWithNormW < currentType>,
		"Householder with normal w in-place",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};*/
	
	methods[3] = {
		new GivensMethodBasic<currentType>(),
		"Givens basic version",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};
	
	/*methods[4] = {
		new GivensMethodQRInOneMatrix<currentType>(),
		"Givens with less memory accesses version",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};
	
	methods[5] = {
		new GivensMethodSIMD<currentType>(),
		"Givens SIMD",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};*/

	for (const auto& method : methods)
	{
		std::cout << method.second.description << std::endl;
		for (auto&& size : method.second.sizes)
		{
			QR_decomposition_test_with_generated_matrix(method.second.solver, size, false);
		}
	}

	std::cout << "all QR tests passed\n";
}

// »спользование: QR разложение матрицы ’ессенберга с использованием вращений
// ћатрицы ’ессенберга:
// * верхн€€ (квадратна€ матрица, у которой все элементы лежащие ниже первой поддиагонали равны нулю, т.е. m[i,j]=0 дл€ любого i>j+1)
// * нижн€€ (при транспонировании получаетс€ верхн€€ матрица ’ессенберга)
void hassenberg_tests()
{
	std::cout << "hassenberg matrix decomposition\n";
	for (auto&& N : { 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 })
	{
		auto res = generateHassenbergMatrix<currentType>(N, true, true);
		// printMatrix(res);
		std::cout << "test size: " << N << std::endl;
		std::vector<std::vector<currentType>> Q(N, std::vector<currentType>(N, 0)), R(N, std::vector<currentType>(N, 0));

		double start = omp_get_wtime();
		IQRSolver<currentType>* solver = new GivensMethodQRInOneMatrix<currentType>();
		solver->QR_decomposition(res, Q, R);
		double end = omp_get_wtime();

		std::cout << "time: " << end - start << std::endl;
		std::cout << "abs error: " << Fnorm(substractMatrix(multiplyMatrix(Q, R), res)) << std::endl;
		std::cout << "rel error: " << Fnorm(substractMatrix(multiplyMatrix(Q, R), res)) / Fnorm(res) << std::endl << std::endl;

		writeMatrixToFile<currentType>("test_data\\matrixA_" + std::to_string(N) + ".txt", res);
		writeMatrixToFile<currentType>("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile<currentType>("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
	std::cout << "all hassenberg tests passed\n\n";
}

int main()
{
	QRtests();
	// hassenberg_tests();
	// std::vector < std::vector <currentType>> A = { {1,2,3}, {4, 5, 6}, {7, 8, 9} }, R, Q;
	//for (auto&& N : { 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 })
	//{
	//	std::vector<std::vector<currentType>> A = generateMatrix<currentType>(N, N);
	//	std::vector<std::vector<currentType>> Q(N, std::vector<currentType>(N, 0)), R(N, std::vector<currentType>(N, 0));
	//	double start = omp_get_wtime();
	//	QR_decomposition(A, Q, R);
	//	double end = omp_get_wtime();

	//	std::cout << "size: " << N << std::endl;
	//	std::cout << "time: " << end - start << std::endl;
	//	std::cout << "abs error: " << Fnorm<currentType>(substractMatrix<currentType>(multiplyMatrix<currentType>(Q, R), A)) << std::endl;
	//	std::cout << "rel error: " << Fnorm<currentType>(substractMatrix<currentType>(multiplyMatrix<currentType>(Q, R), A)) / Fnorm<currentType>(A) << std::endl << std::endl;
	//	//printResultWithExpected(A, Q, R);
	//}
	return 0;
}

