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

	methods[0] = {
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
	};

	methods[3] = {
		new GivensMethodBasic<currentType>(),
		"Givens basic version",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

	methods[4] = {
		new GivensMethodQRInOneMatrix<currentType>(),
		"Givens with less memory accesses version",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

	methods[5] = {
		new GivensMethodSIMD<currentType>(),
		"Givens SIMD",
		{100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

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

template <typename T>
void rotate_rows(std::vector<std::vector<T>>& A, int i, int j, double c, double s)
{
	int N = A.size();
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
	for (auto k = j; k < N; k++)
	{
		// мен€ютс€ только две строки -- i и j
		double temp = A[j][k] * c - A[i][k] * s;
		A[i][k] = A[j][k] * s + A[i][k] * c;
		A[j][k] = temp;
		// T temp = static_cast<T>(R[j][k] * c - R[i][k] * s);
		// R[i][k] = static_cast<T>(R[j][k] * s + R[i][k] * c);
		// R[j][k] = temp;
	}
}

template <typename T>
void rotate_cols(std::vector<std::vector<T>>& A, int i, int j, double c, double s)
{
	int N = A.size();
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
	for (auto k = 0; k < N; k++)
	{
		// мен€ютс€ только два столбца -- i и j
		double temp = c * A[k][j] - s * A[k][i];
		A[k][i] = s * A[k][j] + c * A[k][i];
		A[k][j] = temp;
	}
}

double count_tau(double c, double s)
{
	double tau = s / (1.0 + c);
	return tau;
}

void count_cs(double ajj, double aij, double& c, double& s)
{
	double sqrt_val = std::sqrt(ajj * ajj + aij * aij);
	c = ajj / sqrt_val;
	s = -aij / sqrt_val;
	// double sqrt_val = std::sqrt(Rjj * Rjj + Rij * Rij);
	// double c = Rjj / sqrt_val;
	// double s = -Rij / sqrt_val;
}

void count_c(double tau, double& c, double& s)
{
	tau = tau;
	double tau2 = tau * tau;
	double denom = 1.0 + tau2;
	c = (1.0 - tau2) / denom;
	s = 2.0 * tau / denom;
}

template <typename T>
void qr_decompos_givens(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R)
{
	auto N = A.size();
	Q = std::vector<std::vector<T>>(N, std::vector<T>(N, 0));
	R = A;

	for (int i = 0; i < N; ++i)
		Q[i][i] = T(1);

	for (int j = 0; j < N - 1; ++j) {
		for (int i = j + 1; i < N; ++i) {
			// приведение к double дл€ точных проверок и вычислений
			double Rjj = static_cast<double>(R[j][j]);
			double Rij = static_cast<double>(R[i][j]);

			// проверки малости (пороги дл€ double)
			if (std::abs(Rij) < 1e-11)
				continue;

			double c, s;
			count_cs(Rij, Rjj, c, s);

			// вычисление матрицы R путЄм вращени€ строк
			rotate_rows(R, i, j, c, s);

			// сохранение коэффициентов вращени€
			R[i][j] = static_cast<T>(count_tau(c, s));
		}
	}

	for (int j = 0; j < N - 1; ++j)
	{
		for (int i = j + 1; i < N; ++i)
		{
			if (std::abs(R[i][j]) < 1e-11) continue;
			// восстановление коэффициентов вращени€
			double tau = static_cast<double>(R[i][j]);
			double c, s;
			count_cs(tau, c, s);

			// вычисление матрицы Q путЄм вращени€ столбцов
			rotate_cols(Q, i, j, c, s);
			R[i][j] = T(0);
		}
	}
}

void test_rotations()
{
	std::vector<std::vector<double>> A = { {1., 2., 3.}, {4., 5., 6.}, {7., 8., 10.} };
	auto N = A.size();

	std::cout << "matrix:";
	printMatrix(A);
	std::cout << std::endl;

	{
		std::vector<std::vector<double>> R(A), Q(N, std::vector<double>(N, 0));
		auto solver = new GivensMethodBasic<currentType>();
		solver->QR_decomposition(A, Q, R);

		std::cout << "expec R:";
		printMatrix(R);
		std::cout << std::endl;

		std::cout << "expec Q:";
		printMatrix(Q);
		std::cout << std::endl;
	}

	{
		std::vector<std::vector<double>> R(A), Q(N, std::vector<double>(N, 0));
		for (auto i = 0; i < N; i++)
			Q[i][i] = 1;

		std::cout << "rotations:\n";
		std::vector<std::vector<int>> tasks = { {2, 0}, {1, 0}, { 2, 1 } };//{ {3, 0}, {3, 1}, {2, 0}, {3, 2}, {2, 1}, {1, 0} };
		for (auto&& rows : tasks)
		{
			int& i = rows[0];
			int& j = rows[1];
			std::cout << "(" << i << "," << j << ")\n";

			// проверки малости (пороги дл€ double)
			if (std::abs(R[i][j]) < 1e-11)
				continue;

			double c, s;
			count_cs(R[j][j], R[i][j], c, s);
			rotate_rows(R, i, j, c, s);

			std::cout << "R:";
			printMatrix(R);
			std::cout << std::endl;

			rotate_cols(Q, i, j, c, s);
		}

		std::cout << "R:";
		printMatrix(R);
		std::cout << std::endl;

		std::cout << "Q:";
		printMatrix(Q);
		std::cout << std::endl;

		std::cout << "Q*R:";
		printMatrix(multiplyMatrix<currentType>(Q, R));
		std::cout << std::endl;

		std::cout << "abs error: " << Fnorm<currentType>(substractMatrix<currentType>(multiplyMatrix<currentType>(Q, R), A)) << std::endl;
		std::cout << "rel error: " << Fnorm<currentType>(substractMatrix<currentType>(multiplyMatrix<currentType>(Q, R), A)) / Fnorm<currentType>(A) << std::endl << std::endl;
	}
}


int main()
{
	// QRtests();
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
	test_rotations();
	return 0;
}
