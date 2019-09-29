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
	LU::ComputeLU(m2, L, U);
}