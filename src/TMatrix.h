#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

#include "TVector.h"

template <typename T> class TMatrix {
public:
    TMatrix<T>() : _matrix() {}

    TMatrix<T>(size_t n, size_t m, T def_value = T()) : _matrix(n* m, def_value), _n(n), _m(m) {}

    TMatrix<T>(const std::initializer_list<TVector<T>>& l) : _matrix(l.size() * l.begin()->Size()), _n(l.size()), _m(l.begin()->Size()) {
        size_t index = 0;
        for (auto&& row : l) {
            for (size_t j = 0; j < row.Size(); ++j) {
                _matrix[index] = row[j];
                ++index;
            }
        }
    }

    TMatrix<T>(const TMatrix<T>& other) : _matrix(other._matrix), _n(other._n), _m(other._m) {}

    TMatrix<T> operator-(const TMatrix<T>& m) const {
        TMatrix<T> res(*this);
        auto N = _matrix.Size();
        for (size_t index = 0; index < GetN(); ++index) {
            res._matrix[index] -= m._matrix[index];
        }
        return res;
    }

    TMatrix<T> operator*(T a) const {
        TMatrix<T> res(_n, _m, 0);
        for (size_t index = 0; index < _matrix.Size(); ++index) {
            res.get(index/_m, index%_m) = _matrix[index] * a;
        }
        return res;
    }

    TVector<T> operator*(const TVector<T>& v) const {
        TVector<T> res(_n, 0);
        for (size_t index = 0; index < _matrix.Size(); ++index) {
            size_t i = index / _m, j = index % _m;
            res[i] += _matrix[index] * v[j];
        }
        return res;
    }

    TMatrix<T> operator*(const TMatrix<T>& m) const {
        assert(GetN() == m.GetN());
        auto N = GetN();
        TMatrix<T> res(N, N, 0);
        int s = 50;
        for (int jj = 0; jj < N; jj += s) {
            for (int kk = 0; kk < N; kk += s) {
                for (int i = 0; i < N; i++) {
                    for (int j = jj; j < ((jj + s) > N ? N : (jj + s)); j++) {
                        T temp = 0;
                        for (int k = kk; k < ((kk + s) > N ? N : (kk + s)); k++) {
                            temp += _matrix[i * _n + k] * m.get(k, j);
                        }
                        res.get(i, j) += temp;
                    }
                }
            }
        }
        return res;
    }

    TMatrix<T>& operator=(const TMatrix<T>& other) {
        _matrix = other._matrix;
        _n = other._n;
        _m = other._m;
        return *this;
    }

    bool operator==(const TMatrix<T>& m) const {
        const float eps = 0.000001;
        for (size_t index = 0; index < _matrix.size(); ++index) {
            if (abs(_matrix[index] - m._matrix[index]) > eps) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const TMatrix<T>& v) const {
        return !(*this == v);
    }

    bool IsZero() const {
        if (_n == 0) {
            return true;
        }
        for (size_t index = 0; index < _matrix.Size(); ++index) {
            if (_matrix[index] != 0) {
                return false;
            }
        }
        return true;
    }

    auto GetN() const {
        return _n;
    }

    auto GetM() const {
        return _m;
    }

    T& get(size_t i, size_t j) {
        return _matrix[i * _n + j];
    }

    const T& get(size_t i, size_t j) const {
        return _matrix[i * _n + j];
    }

    friend std::ostream& operator<<(std::ostream& os, const TMatrix<T>& m) {
        os << "[";
        auto N = m._matrix.Size();
        auto M = m._matrix[0].Size();
        for (auto i = 0; i < N; ++i) {
            std::cout << std::left << std::setw(5);
            for (auto j = 0; j < M - 1; ++j) {
                os << m._matrix[i * _n + j] << ", ";
            }
            os << m._matrix[i * _n + M - 1];
            if (i == N - 1) {
                os << "]";
            }
            os << std::endl;
        }
        return os;
    }


private:
    TVector<T> _matrix;
    size_t _n, _m;
};

template <typename T>
TMatrix<T> generate_matrix(int N, int M) {
    TMatrix<T> res(N, M);
    std::srand(std::time(0));
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            res.get(i, j) = std::rand() % 20 - 10;
        }
    }
    return res;
}

template <typename T>
TMatrix<T> operator*(TVector<T> v1, TVector<T> v2) {
    auto N = v1.Size();
    auto M = v2.Size();
    TMatrix<T> res(N, N, 0);
    for (auto i = 0; i < N; ++i) {
        for (auto j = 0; j < M; ++j) {
            res.get(i, j) = v1[i] * v2[j];
        }
    }
    return res;
}
