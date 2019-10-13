#include <cstddef>
#include <fstream>
#include <iostream>
namespace LU{

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

		void operator = (const Matrix &Other);
		double* &operator[](std::size_t index);
		double* operator*(double* vector);

	private:
		std::size_t    dimension_;
		std::ifstream  fin_;
		double **      matrix_;
	};

	void ComputeLU(Matrix & A, Matrix & L, Matrix & U);	
};

void FillVector(double* vector, std::size_t dimention);
void PrintVector(double* vector, std::size_t dimension);
double* FindingY(LU::Matrix & L, double* b);
double* FindingX(LU::Matrix & U, double* y);
double* CalculateDifference(double* vector1, double* vector2, std::size_t dimension);
double CalculateNorm(double* vector, size_t dimension);

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

#include <cmath>
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

double* FindingY(LU::Matrix & L, double* b)
{
	size_t dimension = L.GetDimension();
	double* res = new double[dimension];
		/*for(uint32_t i(0); i < dimension; ++i)
			res[i] = 0.0;*/
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
		/*for(uint32_t i(0); i < dimension; ++i)
			res[i] = 0.0;*/

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

double* CalculateDifference(double* vector1, double* vector2, std::size_t dimension)
{
	for(uint32_t i(0); i < dimension; ++i)
			vector1[i] -= vector2[i];
	return vector1;
}



int main(int argc, char* argv[])
{
	if(argc == 0)
		std::cout << "Input filename!" << std::endl;
	// std::cout << "Hello!" << std::endl;
	LU::Matrix matrix;
	// LU::Matrix m1(2);
	LU::Matrix m2(argv[1]);

	m2.Read();
	// m2.Print();
	std::size_t dimension = m2.GetDimension();

	LU::Matrix L(dimension);
	LU::Matrix U(dimension);
	LU::ComputeLU(m2, L, U);

	double* vectorX = new double[dimension]; // 1 2 3 4 5
	FillVector(vectorX, dimension);
	std::cout << "Vector X = " << std::endl;
	PrintVector(vectorX, dimension);
	
	double* vectorB = m2 * vectorX;
	
	std::cout << "Vector B = " << std::endl;
	PrintVector(vectorB, dimension);

	double* vectorY = FindingY(L, vectorB);
	/*std::cout << "Vector Y = " << std::endl;
	PrintVector(vectorY, dimension);*/

	double* vectorXX = FindingX(U, vectorY);
	std::cout << "Vector Y = " << std::endl;
	PrintVector(vectorXX, dimension);
	// PrintVector(vectorX, dimension);

	double *VectorDif = CalculateDifference(vectorXX, vectorX, dimension);
	// PrintVector(VectorDif, dimension);
	double norm = CalculateNorm(VectorDif, dimension);
	std::cout << "norm(x - y) = " << norm << std::endl;
}