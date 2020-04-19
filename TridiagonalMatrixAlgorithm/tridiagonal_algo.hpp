#ifndef TRIDIAGONAL_ALGO_HPP
#define TRIDIAGONAL_ALGO_HPP

#include <cstddef>
#include <vector>
#include <cmath>

static const double A = 1.0;
static const double B = 2.0 * exp(1);

class Matrix
{
public:
	Matrix(){};
	Matrix(std::size_t dimension);
	~Matrix();
	std::size_t GetDim() const;
	void PrintVector();
	bool TridiagonalAlgo();
private:
	std::size_t dimension_;
	std::vector<double> up_;
	std::vector<double> middle_;
	std::vector<double> down_;
};

bool TridiagonalAlgo(Matrix & matrix);

#endif // TRIDIAGONAL_ALGO_HPP