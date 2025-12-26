#pragma once

#include <vector>

template <typename T>
class IQRSolver {
public:
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) = 0;
};
