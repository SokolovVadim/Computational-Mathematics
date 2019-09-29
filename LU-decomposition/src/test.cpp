#include <iostream>
#include "../include/LU.hpp"


int main()
{
	std::cout << "Hello!" << std::endl;
	LU::Matrix matrix;
	LU::Matrix m1(2);
	LU::Matrix m2("test/mytest.txt");

	m2.Read();
	m2.Print();
	std::size_t dimension = m2.GetDimension();

	LU::Matrix L(dimension);
	LU::Matrix U(dimension);
	LU::ComputeL(m2, L, U);
	/*std::ifstream fin;
	fin.open("test/mytest.txt");
	if(!fin)
		{
			std::cerr << "Unable to open file " << std::endl;
			exit(EXIT_FAILURE); /// STYLE IS NOTHING
		}
		*/
}