#include <iostream>
#include "../include/LU.hpp"

int main(int argc, char* argv[])
{
	if(argc == 0)
		std::cout << "Input filename!" << std::endl;
	std::cout << "Hello!" << std::endl;
	LU::Matrix matrix;
	LU::Matrix m1(2);
	LU::Matrix m2(argv[1]);

	m2.Read();
	m2.Print();
	std::size_t dimension = m2.GetDimension();

	LU::Matrix L(dimension);
	LU::Matrix U(dimension);
	LU::ComputeLU(m2, L, U);

	double* vectorX = new double[dimension];
	FillVector(vectorX, dimension);
	// PrintVector(vectorX, dimension);
	double* vectorB = m2 * vectorX;
	std::cout << "Vector B = " << std::endl;
	PrintVector(vectorB, dimension);

	double* vectorY = FindingY(L, vectorB);
	std::cout << "Vector Y = " << std::endl;
	PrintVector(vectorY, dimension);

	double* vectorZ = FindingX(U, vectorY);
	std::cout << "Vector X = " << std::endl;
	PrintVector(vectorZ, dimension);
}