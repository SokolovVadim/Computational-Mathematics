#include <cstddef>
// #include <iostream>
#include <fstream>
namespace LU{

	class Matrix
	{
	public:
		Matrix();
		Matrix(std::size_t dimention);
		Matrix(const std::string & filename);

		~Matrix();

		void Print();
		bool ReadAndSetDimention();
		bool Read();

	private:
		std::size_t    dimention_;
		std::ifstream  fin_;
		double **      matrix_;
	};
};