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


#include <string>
#include <fstream>
#include <iomanip>


namespace Newton{
	
	const double MAX_NORM = 100000.0;
	const double MIN_NORM = 1e-6;

	void 		Construct_diff_matrix	(LU::Matrix & matrix_a, LU::Matrix & matrix_m, double* vector_exp);
	void        Construct_vector_c      (LU::Matrix & matrix_a, double* vector_exp, double* vector_c, std::size_t dimension);
	void        Linear_system_solution  (LU::Matrix & matrix_a, double* vector_u);
	
};

void 		FillVector 				(double* vector, std::size_t dimension);
void 		PrintVector				(double* vector, std::size_t dimension);
double 		CalculateNorm			(double* vector, size_t dimension);
double* 	CalculateDifference		(double* vector1, double* vector2, std::size_t dimension);
double 		CalculateDistance		(double* vector1, double* vector2, std::size_t dimension);
void 		FillVectorExp			(double* vector_u, double* vector_exp, std::size_t dimension);



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
			vector_c[i] = 0;
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
		vector_u[0] = 1.001;
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
			std::cout << "Solution Mg = c:" << std::endl;
			PrintVector(vectorXX, dimension);

			// LU decomposition ends here

			/*std::cout << "vector u:" << std::endl;
			PrintVector(vector_u, dimension);
			std::cout << "vector g:" << std::endl;
			PrintVector(vector_g, dimension);*/

			for(std::size_t i(0); i < dimension; ++i)
				vector_g[i] = vector_u[i]; // u_n

			for(std::size_t i(0); i < dimension; ++i)
				vector_u[i] += vectorXX[i]; // u_n+1
			std::cout << "vector g:" << std::endl;
			PrintVector(vector_g, dimension);
			
			std::cout << std::endl << "Nearest solution: vector u:" << std::endl;
			std::cout << "****************************" << std::endl;
			PrintVector(vector_u, dimension);
			std::cout << "****************************" << std::endl;
			

			norm = CalculateDistance(vector_g, vector_u, dimension);
			std::cout << "-------------------------------------------------------------------------------" << std::endl;
			std::cout << "step = " << step << " norm = " << norm << std::endl;
			std::cout << "-------------------------------------------------------------------------------" << std::endl;
			step++;



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


#include <iostream>


int main(int argc, char* argv[]	)
{
	std::cout << "Hello!" << std::endl;

	LU::Matrix matrix_a(argv[1]);

	matrix_a.Read();
	matrix_a.Print();

	std::size_t dimension = matrix_a.GetDimension();
	std::cout << "dimension = " << dimension << std::endl;

	

	// construct_diff_matrix()
	// construct vector c
	// cycle: solve equation and calculate norm
	double* vector_u = new double[dimension];
	FillVector(vector_u, dimension);

	/*
	double* vector_c = new double[dimension];
	Newton::Construct_vector_c(matrix_a, vector_u, vector_c, dimension);
	std::cout << "Vector_c:" << std::endl;
	PrintVector(vector_c, dimension);
	*/

	/*
	FillVectorExp(vector_exp, dimension);
	std::cout << "Vector_exp:" << std::endl;
	PrintVector(vector_exp, dimension);
	*/

	// Newton::Matrix matrix_m(dimension);

	
	Newton::Linear_system_solution(matrix_a, vector_u);

	
	delete[] vector_u;
	// delete[] vector_c;
}



