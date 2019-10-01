#include "../include/LU.hpp"
#include <iostream>
#include <sstream>

namespace LU
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

	void Matrix::operator=(const Matrix &Other)
	{
		std::size_t n = this->GetDimension();
		for(uint32_t i(0); i < n; ++i)
			for(uint32_t j(0); j < n; ++j)
				this->matrix_[i][j] = Other.matrix_[i][j];
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

	void ComputeLU(Matrix & A, Matrix & L, Matrix & U)
	{
		
		U = A;
		L.FillWithZeoes();

		/*
		std::cout << "L" << std::endl;
		L.Print();
		std::cout << "U" << std::endl;
		U.Print();
		*/
		std::size_t dimension = A.GetDimension();

		// Prepare L matrix

		for(uint32_t i = 0; i < dimension; i++)
			for(uint32_t j = i; j < dimension; j++)
			{
				if(i == j) // dioganal elements == 1
					L[i][j] = 1;
				else
					L[j][i] = U[j][i] / U[i][i];
			}

		/*
		std::cout << "L" << std::endl;
		L.Print();
		*/

		// Fulfill U and L matrix
	
		for(uint32_t k = 1; k < dimension; k++)
		{
			for(uint32_t i = k - 1; i < dimension; i++)
				for(uint32_t j = i; j < dimension; j++)
				{
					if(i == j)
						L[j][i] = 1;
					else
						L[j][i] = U[j][i] / U[i][i];
				}

			for(uint32_t i = k; i < dimension; i++)
				for(uint32_t j = k - 1; j < dimension; j++)
					U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
		}

		/*std::cout << "L" << std::endl;
		L.Print();
		std::cout << "U" << std::endl;
		U.Print();*/

	}


	// ----------------------------------------------------------------

	

}