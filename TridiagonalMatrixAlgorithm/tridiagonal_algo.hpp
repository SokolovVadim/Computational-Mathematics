#ifndef TRIDIAGONAL_ALGO_HPP
#define TRIDIAGONAL_ALGO_HPP

#include <cstddef>
#include <vector>
#include <cmath>
#include <fstream>

static const double A = 1.0;
static const double B = 2.0 * exp(1);

class Matrix
{
public:
	Matrix(){};
	Matrix(std::size_t dimension, const std::string & filename);
	~Matrix();
	std::size_t GetDim() const;
	void PrintVector();
	void SaveResults(std::vector<double> & approximate_solution);
	std::vector<double> TridiagonalAlgo();
	double CalculateNorm(std::vector<double> & approximate_solution);
private:
	std::size_t dimension_;
	std::vector<double> up_;
	std::vector<double> middle_;
	std::vector<double> down_;
	std::ofstream fout_;
};

bool TridiagonalAlgo(Matrix & matrix);

#endif // TRIDIAGONAL_ALGO_HPP