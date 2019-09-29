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

	private:
		std::size_t    dimension_;
		std::ifstream  fin_;
		double **      matrix_;
	};

	void ComputeU(Matrix & A, Matrix & L, Matrix & U);
	void ComputeL(Matrix & A, Matrix & L, Matrix & U);
};