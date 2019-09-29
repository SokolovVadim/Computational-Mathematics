#include <cstddef>
// #include <iostream>
#include <fstream>
namespace LU{

	class Matrix
	{
	public:
		Matrix();
		Matrix(std::size_t dimention);
		Matrix(std::size_t dimention, const std::string & filename);

		~Matrix();

		void Print();
		bool Read();

	private:
		std::size_t    dimention_;
		std::ifstream  fin_;
		double **      matrix_;
	};
};