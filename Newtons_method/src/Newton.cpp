#include "../include/Newton.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

namespace Newton
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

	void Construct_diff_matrix(Newton::Matrix & matrix_a, Newton::Matrix & matrix_m, double* vector_exp)
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

	void Construct_vector_c(Newton::Matrix & matrix_a, double* vector_exp, double* vector_c, std::size_t dimension)
	{
		// c[i] = -f[i]
		for(std::size_t i(0); i < dimension; ++i)
		{
			for(std::size_t j(0); j < dimension; ++j)
				vector_c[i] -= matrix_a[i][j] * vector_exp[j];
			vector_c[i] += exp(-vector_exp[i]);
		}
	}

	void Linear_system_solution  (Matrix & matrix_a, double* vector_exp, double* vector_c, double* vector_init)
	{
		std::size_t dimension = matrix_a.GetDimension();
		Matrix matrix_m(dimension);
		double norm(MIN_NORM);
		double* vector_g = new double[dimension];

		// std::cout << "min norm: " << std::setprecision(12) << MIN_NORM + 1.0 << std::endl;
		for(std::size_t i(0); i < dimension; ++i)
			vector_g[i] = vector_init[i];
		while(norm > MIN_NORM)
		{

		}

		Newton::Construct_diff_matrix(matrix_a, matrix_m, vector_exp);
		matrix_m.Print();

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
		std::cout << std::setprecision(6) << vector[i] << "\t"; // std::fixed << std::setprecision() << ......
	std::cout << std::endl;
}

double CalculateNorm(double* vector, size_t dimension)
{
	double norm(0.0);
	for(uint32_t i(0); i < dimension; ++i)
		norm += vector[i] * vector[i];
	return sqrt(norm);
}

void FillVectorExp(double* vector, std::size_t dimension)
{
	for(uint32_t i(0); i < dimension; ++i)
	{
		vector[i] = exp(-vector[i]);
	}
	// return vector;
}
