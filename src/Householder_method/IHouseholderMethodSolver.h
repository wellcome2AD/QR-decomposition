#pragma once

#include <vector>

#include "../IQRSolver.h"

template <typename T>
class IHouseholderMethodSolver : public IQRSolver<T> {
public:
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) = 0;
};
