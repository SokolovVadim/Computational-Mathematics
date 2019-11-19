#include "../include/Newton.hpp"
#include <iostream>

int main(int argc, char* argv[]	)
{
	std::cout << "Hello!" << std::endl;

	Newton::Matrix matrix(argv[1]);

	matrix.Read();
	matrix.Print();

	std::size_t dimension = matrix.GetDimension();
	std::cout << "dimension = " << dimension << std::endl;
	
}