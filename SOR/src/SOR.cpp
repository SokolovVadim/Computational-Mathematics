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
				std::cout << matrix_[i][j] << " ";
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


	int Successive_Overrelaxation_algorithm()
	{
		
	}

};



void FillVector(double* vector, std::size_t dimension)
{
	// std::cout << "Vector:" << std::endl;
	for(uint32_t i(0); i < dimension; ++i)
		vector[i] = i + 1;
}

void PrintVector(double* vector, std::size_t dimension)
{
	for(uint32_t i(0); i < dimension; ++i)
		std::cout << vector[i] << std::endl;
}

double CalculateNorm(double* vector, size_t dimension)
{
	double norm(0.0);
	for(uint32_t i(0); i < dimension; ++i)
		norm += vector[i] * vector[i];
	return sqrt(norm);
}