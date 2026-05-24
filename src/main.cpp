#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <map>
#include <vector>
#include <thread>

#include "IQRSolver.h"

#include "HouseholderMethod/HouseholderBasic.h"
#include "HouseholderMethod/HouseholderWithoutMatrixMultiplication.h"
#include "HouseholderMethod/HouseholderWithNormVInplace.h"

#include "GivensMethod/GivensBasic.h"
#include "GivensMethod/GivensQRInOneStruct.h"
#include "GivensMethod/GivensVectorized.h"

#include "Matrix/MatrixOperations.h"

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
void QR_decomposition_test(IQRSolver<T>* solver, const std::vector<std::vector<T>> &A, int N, bool writeToFile)
{
	std::vector<std::vector<T>> Q(N, std::vector<T>(N, 0)), R(N, std::vector<T>(N, 0));

	double start = omp_get_wtime();
	solver->QR_decomposition(A, Q, R);
	double end = omp_get_wtime();

	//std::cout << end - start << std::endl;
	std::cout << "test size: " << N << std::endl;
	std::cout << "time: " << end - start << std::endl;
	std::cout << "abs error: " << Fnorm<T>(substractMatrix<T>(multiplyMatrix<T>(Q, R), A)) << std::endl;
	std::cout << "rel error: " << Fnorm<T>(substractMatrix<T>(multiplyMatrix<T>(Q, R), A)) / Fnorm<T>(A) << std::endl;
	std::cout << std::endl;

	//printResultWithExpected(A, Q, R);

	if (writeToFile) {
		writeMatrixToFile<T>("test_data\\matrixA_" + std::to_string(N) + ".txt", A);
		writeMatrixToFile<T>("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile<T>("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
}

template <typename T>
void QR_decomposition_test_with_randomly_generated_matrix(IQRSolver<T>* solver, int N, bool writeToFile)
{
	std::vector<std::vector<T>> A = generateMatrix<T>(N, N);
	//std::vector<std::vector<T>> A = { {1,2,3},{4,5,6},{7,8,9} };
	QR_decomposition_test<T>(solver, A, N, writeToFile);
}

template <typename T>
void QR_decomposition_test_with_hessenberg_matrix(IQRSolver<T>* solver, int N, bool upper, bool lower, bool writeToFile)
{
	std::vector<std::vector<T>> A = generateHessenbergMatrix<currentType>(N, upper, lower);
	QR_decomposition_test<T>(solver, A, N, writeToFile);
}

struct testParams
{
	IQRSolver<currentType>* solver;
	std::string description;
	std::vector<int> sizes;
};

void QRtests()
{
	std::cout << "QR tests with random matrix\n";

	auto methods = std::map<int, testParams>{};

	//methods[0] = {
	//	new HouseholderMethodBasic<currentType>(),
	//	"Householder basic version",
	//	{ 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000},
	//};

	//methods[1] = {
	//	new HouseholderMethodWithoutMatrixMults<currentType>(),
	//	"Householder without matrix multiplications",
	//	{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	//};

	methods[2] = {
		new HouseholderWithNormVInplace < currentType>,
		"Householder with normal w in-place",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

	//methods[3] = {
	//	new GivensMethodBasic<currentType>(),
	//	"Givens basic version",
	//	{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	//};

	//methods[4] = {
	//	new GivensQRInOneStruct<currentType>(),
	//	"Givens with less memory accesses version",
	//	{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	//};

	//methods[5] = {
	//	new GivensVectorized<currentType>(),
	//	"Givens SIMD",
	//	{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	//};

	for (const auto& method : methods)
	{
		std::cout << method.second.description << std::endl;
		for (auto&& size : method.second.sizes)
		{
			QR_decomposition_test_with_randomly_generated_matrix(method.second.solver, size, false);
		}
	}

	std::cout << "all QR tests passed\n";
}

// »спользование: QR разложение матрицы ’ессенберга с использованием вращений
// ћатрицы ’ессенберга:
// * верхн€€ (квадратна€ матрица, у которой все элементы лежащие ниже первой поддиагонали равны нулю, т.е. m[i,j]=0 дл€ любого i>j+1)
// * нижн€€ (при транспонировании получаетс€ верхн€€ матрица ’ессенберга)
void hessenberg_tests()
{
	auto methods = std::map<int, testParams>{};

	methods[0] = {
		new HouseholderWithNormVInplace<currentType>(),
		"Householder with normal w in-place",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

	methods[1] = {
		new GivensVectorized<currentType>(),
		"Givens SIMD",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

	std::cout << "upper hessenberg matrix decomposition\n";
	for (const auto& method : methods)
	{
		std::cout << method.second.description << std::endl;
		for (auto&& size : method.second.sizes)
		{
			QR_decomposition_test_with_hessenberg_matrix<currentType>(method.second.solver, size, true, false, false);
		}
	}

	std::cout << "lower hessenberg matrix decomposition\n";
	for (const auto& method : methods)
	{
		std::cout << method.second.description << std::endl;
		for (auto&& size : method.second.sizes)
		{
			QR_decomposition_test_with_hessenberg_matrix<currentType>(method.second.solver, size, false, true, true);
		}
	}

	std::cout << "all hessenberg tests passed\n\n";
}

void compare_two_methods_tests()
{
	auto methods = std::map<int, testParams>{};
	methods[0] = {
		new HouseholderWithNormVInplace<currentType>(),
		"Householder",
	};

	methods[1] = {
		new GivensVectorized<currentType>(),
		"Givens",
	};

	int N = 5000;
	std::cout << "test size: " << N << std::endl;
	std::cout << std::endl;
	for (const auto& zeros_ratio : { 0.0, 0.1, 0.2, 0.5, 0.7, 0.75, 0.8, 0.9, 0.95, 1.0 })
	{
		std::cout << "zeros' ratio: " << zeros_ratio << std::endl;
		std::vector<std::vector<currentType>> A = generateMatrixWithLowerZeros<currentType>(N, N, zeros_ratio);
		std::cout << std::endl;
		for (const auto& method : methods)
		{
			std::vector<std::vector<currentType>> Q, R;
			auto& solver = method.second.solver;
			double start = omp_get_wtime();
			solver->QR_decomposition(A, Q, R);
			double end = omp_get_wtime();

			std::cout << method.second.description << " time: " << end - start << std::endl;
		}
		std::cout << "------------------------" << std::endl;
	}
	std::cout << "all compares tests passed\n\n";

	/*
	test size: 5000

zeros' ratio: 0

Householder time: 49.117
Givens time: 63.2722
------------------------
zeros' ratio: 0.1

Householder time: 46.8949
Givens time: 63.2339
------------------------
zeros' ratio: 0.2

Householder time: 46.6289
Givens time: 63.3328
------------------------
zeros' ratio: 0.5

Householder time: 50.2441
Givens time: 62.9841
------------------------
zeros' ratio: 0.7

Householder time: 47.0082
Givens time: 63.5073
------------------------
zeros' ratio: 0.75

Householder time: 47.2866
Givens time: 63.1858
------------------------
zeros' ratio: 0.8

Householder time: 47.1892
Givens time: 63.3863
------------------------
zeros' ratio: 0.9

Householder time: 46.7706
Givens time: 62.9441
------------------------
zeros' ratio: 0.95

Householder time: 46.9988
Givens time: 62.8217
------------------------
zeros' ratio: 1

Householder time: 46.782
Givens time: 0.404332
------------------------
all compares tests passed
	*/
}

template <typename T>
void QR_decomposition1(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R)
{
	int N = A.size();
	int sum = 0;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (std::abs(A[i][j]) < 1e-11)
				sum++;

	IQRSolver<T>* solver;
	float ratio = sum / ((N * N - N) / 2.0);
	if (ratio < 0.9)
	{
		solver = new HouseholderWithNormVInplace<currentType>();
		std::cout << "householder method was choosen to solve" << std::endl;
	}
	else
	{
		solver = new GivensVectorized<currentType>();
		std::cout << "givens method was choosen to solve" << std::endl;
	}
	std::cout << std::endl;
	solver->QR_decomposition(A, Q, R);
}

void choose_method_tests()
{
	int N = 5000;

	std::cout << "test size: " << N << std::endl;
	for (const auto& zeros_ratio : { 0.0, 0.1, 0.2, 0.5, 0.7, 0.75, 0.8, 0.9, 0.95, 1.0 })
	{
		std::cout << "zeros' ratio: " << zeros_ratio << std::endl;
		std::vector<std::vector<currentType>> A = generateMatrixWithLowerZeros<currentType>(N, N, zeros_ratio);
		std::vector<std::vector<currentType>> Q, R;

		double start = omp_get_wtime();
		QR_decomposition1(A, Q, R);
		double end = omp_get_wtime();
		std::cout << "time: " << end - start << std::endl;
	}

	std::cout << "all choosing method tests passed\n\n";
}

int main()
{
	// test multiplyMatrix block size
	if (false)
	{
		int N = 2000;
		auto A = generateMatrix<currentType>(N, N);
		std::cout << "test size: " << N << std::endl;
		std::cout << std::endl;
		float min_time = 1000;
		int block_size_with_min_time;
		for (auto&& size : { 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, N / 2, N / 4, N / 5, N / 10 })
		{
			std::cout << "block size: " << size << std::endl;
			double start = omp_get_wtime();
			multiplyMatrix(A, A, size);
			double end = omp_get_wtime();
			float time = end - start;
			std::cout << "time: " << time << std::endl;
			std::cout << std::endl;
			if (min_time > time)
			{
				min_time = time;
				block_size_with_min_time = size;
			}
		}

		std::cout << std::endl;
		std::cout << "min_time: " << min_time << std::endl;
		std::cout << "block size: " << block_size_with_min_time << std::endl;
	}

	QRtests();
	//hessenberg_tests();
	//compare_two_methods_tests();
	//choose_method_tests();
	return 0;
}
