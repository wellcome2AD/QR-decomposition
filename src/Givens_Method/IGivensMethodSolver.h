#pragma once

#include "../IQRSolver.h"
#include "../TMatrix.h"

template <typename T>
class IGivensMethodSolver : public IQRSolver<T> {
public:
	virtual void QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R) = 0;
};
