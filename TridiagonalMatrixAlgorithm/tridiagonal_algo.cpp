#include "tridiagonal_algo.hpp"
#include <iostream>

Matrix::Matrix(std::size_t dimension):
	dimension_(dimension)
{
	double h = 1.0 / dimension_;

	// fill up

	up_.push_back(0);
	for(std::size_t i(0); i < dimension_ - 2; ++i)
		up_.push_back(1 - h / 2.0);

	// fill middle

	middle_.push_back(1);
	for(std::size_t i(0); i < dimension_ - 2; ++i)
		middle_.push_back(-2);
	middle_.push_back(1);

	// fill down
	for(std::size_t i(0); i < dimension_ - 2; ++i)
		down_.push_back(1 + h / 2.0);
	down_.push_back(0);

}

void Matrix::PrintVector()
{
	for (auto i : up_)
	   std::cout << i << ' ';
	std::cout << std::endl;

	for (auto i : middle_)
	   std::cout << i << ' ';
	std::cout << std::endl;

	for (auto i : down_)
	   std::cout << i << ' ';
	std::cout << std::endl;
}

Matrix::~Matrix()
{
	up_.clear();
	middle_.clear();
	down_.clear();
}

std::size_t Matrix::GetDim() const
{
	return this->dimension_;
}


bool TridiagonalAlgo(Matrix & matrix)
{
	matrix.PrintVector();
	
	std::vector<double> new_up;
	std::vector<double> new_right;
	




	return true;
}