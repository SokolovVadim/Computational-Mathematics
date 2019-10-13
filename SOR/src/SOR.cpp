#include "../include/SOR.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

namespace SOR
{
	Matrix::Matrix():
		dimension_(0),
		matrix_(nullptr)
	{
		// std::cout << "Default matrix" << std::endl;
	}

	// ----------------------------------------------------------------

	Matrix::Matrix(std::size_t dimension):
		dimension_(dimension)
	{
		matrix_ = new double*[dimension_];

		// std::cout << "Cicle allocation" << std::endl;

		for(uint32_t i(0); i < dimension_; ++i)
		{
			matrix_[i] = new double[dimension_];
		}

		// std::cout << "Matrix constructed" << std::endl;
	}

	// ----------------------------------------------------------------

	Matrix::Matrix(const std::string & filename)
	{
		// file openning
		fin_.open(filename);
		if(!fin_)
		{
			std::cerr << "Unable to open file " << filename << std::endl;
			exit(EXIT_FAILURE); /// STYLE IS NOTHING
		}

		this->ReadAndSetDimension();

		// memory allocation
		matrix_ = new double*[dimension_];
		for(uint32_t i(0); i < dimension_; ++i)
		{
			matrix_[i] = new double[dimension_];
		}
		// std::cout << "Matrix constructed" << std::endl;
	}

	// ----------------------------------------------------------------

	Matrix::~Matrix()
	{
		for(uint32_t i(0); i < dimension_; ++i)
		{
			delete[] matrix_[i];
		}
		delete[] matrix_;
		fin_.close();
	}

	// ----------------------------------------------------------------

	void Matrix::Print()
	{
		std::cout << "Matrix: " << std::endl;
		for(uint32_t i(0); i < dimension_; ++i)
		{
			for(uint32_t j(0); j < dimension_; ++j)
			{
				std::cout << "\t" << matrix_[i][j];
			}
			std::cout << std::endl;
		}
	}

	// ----------------------------------------------------------------

	bool Matrix::Read()
	{
		std::string str((std::istreambuf_iterator<char>(fin_)),
                 std::istreambuf_iterator<char>());
		// std::cout << str;
		std::stringstream stream(str);
		for(uint32_t i(0); i < dimension_; ++i)
		{
			for(uint32_t j(0); j < dimension_; ++j)
			{
				stream >> matrix_[i][j];
			}
		}
		return true;
	}

	// ----------------------------------------------------------------

	bool Matrix::ReadAndSetDimension()
	{
		size_t dimension(0);
		fin_ >> dimension;
		// std::cout << "dimension read: " << dimension << std::endl;
		this->dimension_ = dimension;
		return true;
	}

	// ----------------------------------------------------------------

	const std::size_t Matrix::GetDimension() const
	{
		return this->dimension_;
	}
	
	// ----------------------------------------------------------------

	void Matrix::FillWithZeoes()
	{
		for(uint32_t i(0); i < dimension_; ++i)
			for(uint32_t j(0); j < dimension_; ++j)
				this->matrix_[i][j] = 0.0;
	}

	// ----------------------------------------------------------------

	double * Matrix::operator*(double* vector)
	{
		double* res = new double[dimension_];
		/*for(uint32_t i(0); i < dimension_; ++i)
			res[i] = 0.0;*/
		for(uint32_t i(0); i < dimension_; ++i)
			for(uint32_t j(0); j < dimension_; ++j)
				res[i] += matrix_[i][j] * vector[j];
		// PrintVector(res, dimension_);
		/*for(uint32_t i(0); i < dimension_; ++i)
			vector[i] = res[i];*/
		return res;
	}

	// ----------------------------------------------------------------

	double* &Matrix::operator[](std::size_t index)
	{
		if(index >= dimension_)
		{
			std::cout << "index is out of range" << std::endl;
		}
		return this->matrix_[index];
	}

	// ----------------------------------------------------------------

	int Successive_Overrelaxation_algorithm( Matrix & matrix, double* vectorB, double* vectorInit,
												double min_residual, double relaxation_parameter)
	{
		std::size_t dimension = matrix.GetDimension();
		double residual = CalculateNorm(CalculateDifference(matrix * vectorInit, vectorB, dimension), dimension); // INIT PROPERLY HERE!!!!!
		std::cout << "Initial residual: " << residual << std::endl;

		uint32_t step(0);
		while(residual > min_residual)
		{
			for(uint32_t i(0); i < dimension; ++i)
			{
				double sigma(0.0);
				for(uint32_t j(0); j < dimension; ++j)
				{
					if(i != j)
						sigma += matrix[i][j] * vectorInit[j];
				}
				vectorInit[i] = (1 - relaxation_parameter) * vectorInit[i] + (relaxation_parameter / matrix[i][i]) * (vectorB[i] - sigma);
			}
			residual = CalculateNorm(CalculateDifference(matrix * vectorInit, vectorB, dimension), dimension);
			std::cout << "Step: " << step << ",\tresidual: " << residual << "\t";
			PrintVector(vectorInit, dimension);
			step++;
		}


		return 0;
	}

	/*
	for j in range(A.shape[1]):
        if j != i:
          sigma += A[i][j] * phi[j]
      phi[i] = (1 - omega) * phi[i] + (omega / A[i][i]) * (b[i] - sigma)
	*/

};


double* CalculateDifference(double* vector1, double* vector2, std::size_t dimension)
{
	for(uint32_t i(0); i < dimension; ++i)
			vector1[i] -= vector2[i];
	return vector1;
}

void FillVector(double* vector, std::size_t dimension)
{
	// std::cout << "Vector:" << std::endl;
	for(uint32_t i(0); i < dimension; ++i)
		vector[i] = i + 1;
}

void PrintVector(double* vector, std::size_t dimension)
{
	for(uint32_t i(0); i < dimension; ++i)
		std::cout << std::fixed << std::setprecision(12) << vector[i] << "\t";
	std::cout << std::endl;
}

double CalculateNorm(double* vector, size_t dimension)
{
	double norm(0.0);
	for(uint32_t i(0); i < dimension; ++i)
		norm += vector[i] * vector[i];
	return sqrt(norm);
}