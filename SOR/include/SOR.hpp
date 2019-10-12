#include <string>
#include <fstream>
namespace SOR{

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
		// double* &operator[](std::size_t index);
		// double* operator*(double* vector);

	private:
		std::size_t    dimension_;
		std::ifstream  fin_;
		double **      matrix_;
	};

	int Successive_Overrelaxation_algorithm();
	
};