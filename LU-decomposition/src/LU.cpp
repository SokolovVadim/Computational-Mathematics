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

	Matrix::Matrix(const std::string & filename)
	{
		// file openning
		fin_.open(filename);
		if(!fin_)
		{
			std::cerr << "Unable to open file " << filename << std::endl;
			exit(EXIT_FAILURE); /// STYLE IS NOTHING
		}

		this->ReadAndSetDimention();

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

	bool Matrix::ReadAndSetDimention()
	{
		size_t dimention(0);
		fin_ >> dimention;
		std::cout << "Dimention read: " << dimention << std::endl;
		this->dimention_ = dimention;
		return true;
	}

	// ----------------------------------------------------------------

}