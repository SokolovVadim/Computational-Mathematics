#include <string>
#include <fstream>
#include <iomanip>
namespace SOR{

	/*const double OMEGA = 1.5; // relaxation parameter
	const double RESIDUAL = 1e-6; // norm of residual
	const double min_residual = 100.0;*/
	const double MAX_NORM = 100000.0;

	class Matrix
	{
	public:
		Matrix();
		Matrix(std::size_t dimension);
		Matrix(const std::string & filename);

		~Matrix();

		void Print();
		bool ReadAndSetDimension();
		bool Read();
		const std::size_t GetDimension() const;
		void FillWithZeoes();

		// void operator = (const Matrix &Other);
		double* &operator[](std::size_t index);
		double* operator*(double* vector);

	private:
		std::size_t    dimension_;
		std::ifstream  fin_;
		double **      matrix_;
	};

	int Successive_Overrelaxation_algorithm(Matrix & matrix, double* vectorB, double* vectorInit,
												double min_residual, double relaxation_parameter);
	int Gauss_Seidel_algorithm( Matrix & matrix, double* vectorB,
								double* vectorInit, double min_norm);
	
};

void FillVector(double* vector, std::size_t dimension);
void PrintVector(double* vector, std::size_t dimension);
double CalculateNorm(double* vector, size_t dimension);
double* CalculateDifference(double* vector1, double* vector2, std::size_t dimension);
double CalculateDistance(double* vector1, double* vector2, std::size_t dimension);


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

	int Gauss_Seidel_algorithm( Matrix & matrix, double* vectorB,
								double* vectorInit, double min_norm)
	{
		std::size_t dimension = matrix.GetDimension();
		double* vectorPrev = new double[dimension];
		double norm = MAX_NORM;
		std::size_t step(0);
		while(norm > min_norm)
		{
			for(uint32_t i(0); i < dimension; ++i)
				vectorPrev[i] = vectorInit[i];
			for(uint32_t i(0); i < dimension; ++i)
			{
				double sum(0.0);
				for(uint32_t j(0); j < i; ++j)
					sum += matrix[i][j] * vectorInit[j];
				for(uint32_t j(i + 1); j < dimension; ++j)
					sum += matrix[i][j] * vectorPrev[j];
				vectorInit[i] = (vectorB[i] - sum) / matrix[i][i];
			}
			norm = CalculateDistance(vectorInit, vectorPrev, dimension);
			std::cout << "step = " << step << ", norm = " << norm << std::endl;
			step++;
		}
		delete[] vectorPrev;
		return 0;
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

int main(int argc, char* argv[]	)
{
	std::cout << "Hello!" << std::endl;

	SOR::Matrix matrix(argv[1]);

	matrix.Read();
	matrix.Print();

	std::size_t dimension = matrix.GetDimension();

	double* vectorB = new double[dimension]; // 1 2 3 4 5

	FillVector(vectorB, dimension);
	std::cout << "VectorB: " << std::endl;
	PrintVector(vectorB, dimension);

	double* vectorInit = new double[dimension]; // 1 2 3 4 5
	
	std::cout << "VectorInit: " << std::endl;
	PrintVector(vectorInit, dimension);

	std::cout << "*************** Successive_Overrelaxation_algorithm ***************" << std::endl;
	
	SOR::Successive_Overrelaxation_algorithm(matrix, vectorB, vectorInit, 1e-12, 1.5);
	std::cout << "Result: " << std::endl;
	PrintVector(vectorInit, dimension);

	// clear initial guess vector
	for(uint32_t i(0); i < dimension; ++i)
	{
		vectorInit[i] = 0.0;
	}

	std::cout << "*************** Gauss_Seidel_algorithm ***************" << std::endl;

	SOR::Gauss_Seidel_algorithm(matrix, vectorB, vectorInit, 1e-12);
	std::cout << "Result: " << std::endl;
	PrintVector(vectorInit, dimension);

	delete[] vectorB;
	delete[] vectorInit;
}