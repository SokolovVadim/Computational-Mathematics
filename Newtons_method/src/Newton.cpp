#include "../include/Newton.hpp"
// #include "../../LU-decomposition/include/LU.hpp"


#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

namespace Newton
{
	void Construct_diff_matrix(LU::Matrix & matrix_a, LU::Matrix & matrix_m, double* vector_exp)
	{
		std::size_t dimension = matrix_a.GetDimension();
		for(std::size_t i(0); i < dimension; ++i)
			for(std::size_t j(0); j < dimension; ++j)
			{
				if(i == j)
					matrix_m[i][j] = matrix_a[i][j] + vector_exp[i];
				else
					matrix_m[i][j] = matrix_a[i][j];
			}
	}

	// ----------------------------------------------------------------

	void Construct_vector_c(LU::Matrix & matrix_a, double* vector_u, double* vector_c, std::size_t dimension)
	{
		// c[i] = -f[i]
		for(std::size_t i(0); i < dimension; ++i)
		{
			for(std::size_t j(0); j < dimension; ++j)
				vector_c[i] -= matrix_a[i][j] * vector_u[j];
			vector_c[i] += exp(-vector_u[i]);
		}
	}

	/*
	vector_u -> e(-vector_u)


	*/

	void Linear_system_solution(LU::Matrix & matrix_a, double* vector_u)
	{
		std::size_t dimension = matrix_a.GetDimension();
		LU::Matrix matrix_m(dimension);
		double norm(MAX_NORM);
		double* vector_c = new double[dimension];
		double* vector_g = new double[dimension];
		double* vector_exp = new double[dimension];

		// std::cout << "min norm: " << std::setprecision(12) << MIN_NORM + 1.0 << std::endl;

		uint32_t step(0);
		
		while(norm > MIN_NORM)
		{
			Newton::Construct_vector_c(matrix_a, vector_u, vector_c, dimension);
			std::cout << "Vector_c:" << std::endl;
			PrintVector(vector_c, dimension);
			
			FillVectorExp(vector_u, vector_exp, dimension);
			std::cout << "Vector_exp:" << std::endl;
			PrintVector(vector_exp, dimension);
			
			/*for(std::size_t i(0); i < dimension; ++i)
				vector_g[i] = vector_init[i];*/
			
			Newton::Construct_diff_matrix(matrix_a, matrix_m, vector_exp);
			matrix_m.Print();

			// LU deconposition starts here
			LU::Matrix L(dimension);
			LU::Matrix U(dimension);
			LU::ComputeLU((matrix_m), L, U); // issue here!!!!!!!!!

			double* vectorY = FindingY(L, vector_c);
	

			double* vectorXX = FindingX(U, vectorY);
			std::cout << "Solution = " << std::endl;
			PrintVector(vectorXX, dimension);

			std::cout << "vector u:" << std::endl;
			PrintVector(vector_u, dimension);
			std::cout << "vector g:" << std::endl;
			PrintVector(vector_g, dimension);

			for(std::size_t i(0); i < dimension; ++i)
				vector_g[i] = vectorXX[i];

			for(std::size_t i(0); i < dimension; ++i)
				vector_u[i] += vector_g[i];
			std::cout << "vector u:" << std::endl;
			PrintVector(vector_u, dimension);
			std::cout << "vector g:" << std::endl;
			PrintVector(vector_g, dimension);

			norm = CalculateDistance(vector_g, vector_u, dimension);
			std::cout << "-------------------------------------------------------------------------------" << std::endl;
			std::cout << "step = " << step << " norm = " << norm << std::endl;
			std::cout << "-------------------------------------------------------------------------------" << std::endl;
			step++;

			// LU decomposition ends here


			/*vector_init = solve_eq(...);
			vector_g += vector_init;*/


			// gonna get vector_g => vector_init += vector_g
			// vector_exp = Fill_vector_exp(vector init)

			// break; // just to debug
		}
		delete[] vector_exp;
		delete[] vector_g;
		delete[] vector_c;
	}

};

double CalculateDistance(double* vector1, double* vector2, std::size_t dimension)
{
	double norm(0.0);
	for(uint32_t i(0); i < dimension; ++i)
	{
		norm = pow(vector1[i] - vector2[i], 2);
	}
	return sqrt(norm);
}

void FillVectorExp(double* vector_u, double* vector_exp, std::size_t dimension)
{
	for(uint32_t i(0); i < dimension; ++i)
	{
		vector_exp[i] = exp(-vector_u[i]);
	}
}

