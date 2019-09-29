#include "../include/LU.hpp"
#include <iostream>


namespace LU
{
	Matrix::Matrix():
		dimention_(0),
		matrix_(nullptr)
	{
		std::cout << "Default matrix" << std::endl;
	}

	// ----------------------------------------------------------------

	Matrix::Matrix(std::size_t dimention):
		dimention_(dimention)
	{
		matrix_ = new double*[dimention_];

		std::cout << "Cicle allocation" << std::endl;

		for(uint32_t i(0); i < dimention_; ++i)
		{
			matrix_[i] = new double[dimention_];
		}

		std::cout << "Matrix constructed" << std::endl;
	}

	// ----------------------------------------------------------------

	Matrix::Matrix(std::size_t dimention, const std::string & filename):
		dimention_(dimention)
	{
		// file openning
		fin_.open(filename);
		if(!fin_)
		{
			std::cerr << "Unable to open file " << filename << std::endl;
			exit(EXIT_FAILURE); /// STYLE IS NOTHING
		}
		// memory allocation
		matrix_ = new double*[dimention_];
		for(uint32_t i(0); i < dimention_; ++i)
		{
			matrix_[i] = new double[dimention_];
		}
		std::cout << "Matrix constructed" << std::endl;
	}

	// ----------------------------------------------------------------

	Matrix::~Matrix()
	{
		for(uint32_t i(0); i < dimention_; ++i)
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
		for(uint32_t i(0); i < dimention_; ++i)
		{
			for(uint32_t j(0); j < dimention_; ++j)
			{
				std::cout << matrix_[i][j];
			}
			std::cout << std::endl;
		}
	}

	// ----------------------------------------------------------------

	bool Matrix::Read()
	{
		return true;
	}

	// ----------------------------------------------------------------

}