#include "tridiagonal_algo.hpp"
#include <iostream>
#include <cmath>

Matrix::Matrix(std::size_t dimension):
	dimension_(dimension)
{
	double h = 1.0 / dimension_;

	// fill up

	up_.push_back(0);
	for(std::size_t i(0); i < dimension_ - 1; ++i)
		up_.push_back(1 - h / 2.0);

	// fill middle

	middle_.push_back(1);
	for(std::size_t i(0); i < dimension_ - 1; ++i)
		middle_.push_back(-2);
	middle_.push_back(1);

	// fill down
	for(std::size_t i(0); i < dimension_ - 1; ++i)
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

double function(double x)
{
	return exp(x);
}

std::vector<double> CalculateRightSide(std::size_t dim)
{
	std::vector<double> right;
	std::cout << "right\n";
	for(std::size_t i(0); i < dim + 1; ++i)
		right.push_back(exp(double(i) / dim));

	return right;
}


bool Matrix::TridiagonalAlgo()
{
	this->PrintVector();
	
	
	std::vector<double> right = CalculateRightSide(this->GetDim());
	
	for (auto i : right)
	   std::cout << i << ' ';
	std::cout << std::endl;

	std::vector<double> new_up;
	std::vector<double> new_right;


	return true;
}