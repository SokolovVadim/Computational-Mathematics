#include "../include/LU.hpp"
#include <iostream>
#include <sstream>

namespace LU
{
	Matrix::Matrix():
		dimension_(0),
		matrix_(nullptr)
	{
		std::cout << "Default matrix" << std::endl;
	}

	// ----------------------------------------------------------------

	Matrix::Matrix(std::size_t dimension):
		dimension_(dimension)
	{
		matrix_ = new double*[dimension_];

		std::cout << "Cicle allocation" << std::endl;

		for(uint32_t i(0); i < dimension_; ++i)
		{
			matrix_[i] = new double[dimension_];
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

		this->ReadAndSetDimension();

		// memory allocation
		matrix_ = new double*[dimension_];
		for(uint32_t i(0); i < dimension_; ++i)
		{
			matrix_[i] = new double[dimension_];
		}
		std::cout << "Matrix constructed" << std::endl;
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
		std::cout << "dimension read: " << dimension << std::endl;
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

	void ComputeL(Matrix & A, Matrix & L, Matrix & U)
	{
		// std::size_t n = A.GetDimension();
		//n = 0;
		U = A;
		L.FillWithZeoes();

		std::cout << "L" << std::endl;
		L.Print();
		std::cout << "U" << std::endl;
		U.Print();

		std::size_t dimension_ = A.GetDimension();

		for(uint32_t i = 0; i < dimension_; i++)
			for(uint32_t j = i; j < dimension_; j++)
				L[j][i] = U[j][i] / U[i][i];

		std::cout << "L" << std::endl;
		L.Print();
	
		for(uint32_t k = 1; k < dimension_; k++)
		{
			for(uint32_t i = k-1; i < dimension_; i++)
				for(uint32_t j = i; j < dimension_; j++)
					L[j][i]=U[j][i]/U[i][i];

			for(uint32_t i = k; i < dimension_; i++)
				for(uint32_t j = k-1; j < dimension_; j++)
					U[i][j]=U[i][j]-L[i][k-1]*U[k-1][j];
		}

		std::cout << "L" << std::endl;
		L.Print();
		std::cout << "U" << std::endl;
		U.Print();

	}


	// ----------------------------------------------------------------

	void ComputeU(Matrix & A, Matrix & L, Matrix & U)
	{

	}


	// ----------------------------------------------------------------

	

}