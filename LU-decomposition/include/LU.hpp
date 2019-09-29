#include <cstddef>
// #include <iostream>
#include <fstream>
namespace LU{

	class Matrix
	{
	public:
		Matrix();
		Matrix(std::size_t dimention);
		~Matrix();

		void Print();

	private:
		std::size_t    dimention_;
		double **      matrix_;
		std::ifstream fin;
	};

	// open file by filename
	size_t ReadDimention();
	// read matrix in Matrix structure

};