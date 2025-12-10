#pragma once

#include "TMatrix.h"

template <typename T>
class IQRSolver {
public:
	virtual void QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R) = 0;
};
