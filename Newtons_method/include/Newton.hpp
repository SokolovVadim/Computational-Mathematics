#include <string>
#include <fstream>
#include <iomanip>
namespace Newton{
	
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

	void 		Construct_diff_matrix	(Matrix & matrix_a, Newton::Matrix & matrix_m, double* vector_exp);
	void        Construct_vector_c      (Matrix & matrix_a, double* vector_exp, double* vector_c, std::size_t dimension);
	void        Linear_system_solution  (Matrix & matrix_a, double* vector_exp, double* vector_c);
	
};

void 		FillVector 				(double* vector, std::size_t dimension);
void 		PrintVector				(double* vector, std::size_t dimension);
double 		CalculateNorm			(double* vector, size_t dimension);
double* 	CalculateDifference		(double* vector1, double* vector2, std::size_t dimension);
double 		CalculateDistance		(double* vector1, double* vector2, std::size_t dimension);
void 		FillVectorExp			(double* vector, std::size_t dimension);
