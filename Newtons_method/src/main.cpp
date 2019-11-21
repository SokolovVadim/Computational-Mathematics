#include "../include/Newton.hpp"
#include <iostream>
#include <vector>

int main(int argc, char* argv[]	)
{
	std::cout << "Hello!" << std::endl;

	Newton::Matrix matrix_a(argv[1]);

	matrix_a.Read();
	matrix_a.Print();

	std::size_t dimension = matrix_a.GetDimension();
	std::cout << "dimension = " << dimension << std::endl;

	std::vector<double> v_right_side(dimension);
	for(auto i: v_right_side)
		std::cout << "v_right_side[i] = "<< i << std::endl;

	// construct_diff_matrix()
	// construct vector c
	// cycle: solve equation and calculate norm
	double* vector_exp = new double[dimension];
	FillVector(vector_exp, dimension);

	double* vector_c = new double[dimension];
	Newton::Construct_vector_c(matrix_a, vector_exp, vector_c, dimension);
	std::cout << "Vector_c:" << std::endl;
	PrintVector(vector_c, dimension);

	FillVectorExp(vector_exp, dimension);
	std::cout << "Vector_exp:" << std::endl;
	PrintVector(vector_exp, dimension);

	// Newton::Matrix matrix_m(dimension);


	Newton::Linear_system_solution(matrix_a, vector_exp, vector_c);

	

	delete[] vector_c;
}