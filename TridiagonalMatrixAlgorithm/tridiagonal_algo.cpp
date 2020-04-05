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
	
	
	std::vector<double> right = CalculateRightSide(dimension_);
	
	for (auto i : right)
	   std::cout << i << ' ';
	std::cout << std::endl;

	std::vector<double> new_down;
	std::vector<double> new_right;
	std::vector<double> y;

	new_down.push_back(down_.front() / middle_.front());
	std::cout << new_down.back() << std::endl;

	new_right.push_back(right.front() / middle_.front());
	std::cout << new_right.back() << std::endl;

	std::cout << "Fill new_down\n";
	for(std::size_t i(1); i < dimension_; ++i)
	{
		new_down.push_back(middle_.at(i) - up_.at(i) * new_down.at(i - 1));
		std::cout << "c[i]' = " << new_down.at(i) << std::endl;
		new_right.push_back((right.at(i) - up_.at(i) * new_right.at(i - 1)) / (middle_.at(i) - up_.at(i) * new_down.at(i - 1)));
		std::cout << "d[i]' = " << new_right.at(i) << std::endl;
	}

	new_right.push_back((right.at(dimension_ - 1) - up_.at(dimension_ - 1) * new_right.at(dimension_ - 2)) / (middle_.at(dimension_ - 1) - up_.at(dimension_ - 1) * new_down.at(dimension_ - 2)));

	y.resize(dimension_);
	y.at(dimension_ - 1) = new_right.at(dimension_ - 1);

	for(int i(dimension_ - 2); i >= 0; --i)
	{
		y.at(i) = (right.at(i) - down_.at(i) * (y.at(i + 1))) / middle_.at(i);
	}

	std::cout << "y\n";
	for (auto i : y)
	   std::cout << i << ' ';
	std::cout << std::endl;


	return true;
}