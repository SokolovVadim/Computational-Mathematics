#include <string>
#include <fstream>
#include <iomanip>
namespace SOR{

	/*const double OMEGA = 1.5; // relaxation parameter
	const double RESIDUAL = 1e-6; // norm of residual
	const double min_residual = 100.0;*/

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
	
};

void FillVector(double* vector, std::size_t dimension);
void PrintVector(double* vector, std::size_t dimension);
double CalculateNorm(double* vector, size_t dimension);
double* CalculateDifference(double* vector1, double* vector2, std::size_t dimension);