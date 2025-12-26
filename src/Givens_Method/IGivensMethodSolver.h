#pragma once

#include <vector>

template <typename T>
class IGivensMethodSolver : public IQRSolver<T> {
public:
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) = 0;
};
