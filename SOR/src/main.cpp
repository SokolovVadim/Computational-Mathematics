#include "../include/SOR.hpp"
#include <iostream>

int main(int argc, char* argv[]	)
{
	std::cout << "Hello!" << std::endl;

	SOR::Matrix matrix(argv[1]);

	matrix.Read();
	matrix.Print();

	std::size_t dimension = matrix.GetDimension();

	double* vectorB = new double[dimension]; // 1 2 3 4 5

	FillVector(vectorB, dimension);
	std::cout << "VectorB: " << std::endl;
	PrintVector(vectorB, dimension);

	double* vectorInit = new double[dimension]; // 1 2 3 4 5
	
	std::cout << "VectorInit: " << std::endl;
	PrintVector(vectorInit, dimension);

	SOR::Successive_Overrelaxation_algorithm(matrix, vectorB, vectorInit, 1e-12, 1.5);
	std::cout << "Result: " << std::endl;
	PrintVector(vectorInit, dimension);

	delete[] vectorB;
	delete[] vectorInit;
}