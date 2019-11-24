#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <string>
#include <fstream>
#include <iomanip>
#include "../../LU-decomposition/include/LU.hpp"

namespace Newton{
	
	const double MAX_NORM = 100000.0;
	const double MIN_NORM = 1e-6;

	void 		Construct_diff_matrix	(LU::Matrix & matrix_a, LU::Matrix & matrix_m, double* vector_exp);
	void        Construct_vector_c      (LU::Matrix & matrix_a, double* vector_exp, double* vector_c, std::size_t dimension);
	void        Linear_system_solution  (LU::Matrix & matrix_a, double* vector_u);
	
};

void 		FillVector 				(double* vector, std::size_t dimension);
void 		PrintVector				(double* vector, std::size_t dimension);
double 		CalculateNorm			(double* vector, size_t dimension);
double* 	CalculateDifference		(double* vector1, double* vector2, std::size_t dimension);
double 		CalculateDistance		(double* vector1, double* vector2, std::size_t dimension);
void 		FillVectorExp			(double* vector_u, double* vector_exp, std::size_t dimension);

#endif // NEWTON_HPP
