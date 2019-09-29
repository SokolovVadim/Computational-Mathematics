#include <iostream>
#include <cmath>
#include "../include/LU.hpp"

void FillVector(double* vector, std::size_t dimension)
{
	std::cout << "Vector:" << std::endl;
	for(uint32_t i(0); i < dimension; ++i)
		vector[i] = i + 1;
}

void PrintVector(double* vector, std::size_t dimension)
{
	for(uint32_t i(0); i < dimension; ++i)
		std::cout << vector[i] << std::endl;
}

double* FindingY(LU::Matrix & L, double* b)
{
	size_t dimension = L.GetDimension();
	double* res = new double[dimension];
		for(uint32_t i(0); i < dimension; ++i)
			res[i] = 0.0;
	for(uint32_t i(0); i < dimension; ++i)
	{
		double sum = 0;
		for(uint32_t j(0); j < i; ++j)
			sum += L[i][j] * res[j];
		res[i] = (b[i] - sum) / L[i][i];
	}
	return res;
}


double* FindingX(LU::Matrix & U, double* y)
{
	size_t dimension = U.GetDimension();
	// std::cout << "dimU = " << dimension << std::endl;
	double* res = new double[dimension];
		for(uint32_t i(0); i < dimension; ++i)
			res[i] = 0.0;

	// std::cout << "I'm here" << std::endl; 
	for(int i = (dimension - 1); i >= 0; i--)
	{
		// std::cout << " i = " << i;
		// std::cout << "I'm here" << std::endl;
		double sum = 0;
		for(int j(dimension - 1); j > i; j--)
			sum += U[i][j] * res[j];
		res[i] = (y[i] - sum) / U[i][i];
		// std::cout << " res[i] = " << res[i];
	}
	return res;
}

double CalculateNorm(double* vector, size_t dimension)
{
	double norm(0.0);
	for(uint32_t i(0); i < dimension; ++i)
		norm += vector[i] * vector[i];
	return sqrt(norm);
}