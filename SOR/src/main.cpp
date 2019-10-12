#include "../include/SOR.hpp"
#include <iostream>

int main(int argc, char* argv[]	)
{
	std::cout << "Hello!" << std::endl;

	SOR::Matrix m(argv[1]);

	m.Read();
	m.Print();
	// std::size_t dimension = m.GetDimension();


}