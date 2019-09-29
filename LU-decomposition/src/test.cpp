#include <iostream>
#include "../include/LU.hpp"


int main()
{
	std::cout << "Hello!" << std::endl;
	LU::Matrix matrix;
	LU::Matrix m1(2);
	LU::Matrix m2(2, "test/mytest.txt");
	/*std::ifstream fin;
	fin.open("test/mytest.txt");
	if(!fin)
		{
			std::cerr << "Unable to open file " << std::endl;
			exit(EXIT_FAILURE); /// STYLE IS NOTHING
		}
		*/
}