#include <iostream>
#include "../include/LU.hpp"

int main(int argc, char* argv[])
{
	if(argc == 0)
		std::cout << "Input filename!" << std::endl;
	// std::cout << "Hello!" << std::endl;
	LU::Matrix matrix;
	// LU::Matrix m1(2);
	LU::Matrix m2(argv[1]);

	m2.Read();
	// m2.Print();
	std::size_t dimension = m2.GetDimension();

	LU::Matrix L(dimension);
	LU::Matrix U(dimension);
	LU::ComputeLU(m2, L, U);

	double* vectorX = new double[dimension]; // 1 2 3 4 5
	FillVector(vectorX, dimension);
	std::cout << "Vector X = " << std::endl;
	PrintVector(vectorX, dimension);
	
	double* vectorB = m2 * vectorX;
	
	std::cout << "Vector B = " << std::endl;
	PrintVector(vectorB, dimension);

	double* vectorY = FindingY(L, vectorB);
	/*std::cout << "Vector Y = " << std::endl;
	PrintVector(vectorY, dimension);*/

	double* vectorXX = FindingX(U, vectorY);
	std::cout << "Vector Y = " << std::endl;
	PrintVector(vectorXX, dimension);
	// PrintVector(vectorX, dimension);

	double *VectorDif = CalculateDifference(vectorXX, vectorX, dimension);
	// PrintVector(VectorDif, dimension);
	double norm = CalculateNorm(VectorDif, dimension);
	std::cout << "norm(x - y) = " << norm << std::endl;
}