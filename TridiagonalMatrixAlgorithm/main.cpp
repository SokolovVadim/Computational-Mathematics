#include "tridiagonal_algo.hpp"
#include <iostream>
#include <climits>

long int ReadArg(char * str);

int main(int argc, char** argv)
{
	if(argc < 2){
		std::cout << "Number of grid intervals required!\n";
		exit(EXIT_FAILURE);
	}
	std::size_t n = ReadArg(argv[1]);
	
	Matrix matrix(n);
	TridiagonalAlgo(matrix);

	return 0;
}

long int ReadArg(char * str)
{
	char* endptr;
	errno = 0;

	long int number = strtol(str, &endptr, 10);

	
	if ((errno == ERANGE && (number == LONG_MAX || number == LONG_MIN)) || (errno != 0 && number == 0)) 
	{
       		perror("strtol");
        	exit(EXIT_FAILURE);
   	}

	if (endptr == str)
	{
        	fprintf(stderr, "Error!\n");
        	exit(EXIT_FAILURE);
   	}
	if (*endptr != '\0')
	{
        	fprintf(stderr, "Error!\n");
        	exit(EXIT_FAILURE);
   	}

	return number;
}