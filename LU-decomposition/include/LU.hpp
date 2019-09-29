#include <cstddef>
// #include <iostream>
#include <fstream>
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
