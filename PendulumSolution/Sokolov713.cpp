#include <iostream>
#include <cmath>
#include <fstream>
#define _USE_MATH_DEFINES

std::ofstream fout("output.txt");

enum INIT_VALUES{
	U0 = 0,
	V0 = 1,
	G = 10,
	LAMBDA = 1,
	L = 1,
	N = 100,
	PERIOD_NUMBER = 3
};

// counting second derivative
double second_derivative(double tau, double coord, double velocity)
{
	double der = (-G * sin(coord) - 2 * LAMBDA * velocity) / (L);
	// std::cout << "second_der = " << der << "\t";
	return der;
}


// 
void EulerSolution(double tau, double * coord, double * velocity){
	coord[0]    = U0;
	velocity[0] = V0;
	for(int i(0); i < PERIOD_NUMBER * N; ++i)
	{
		velocity[i + 1] = velocity[i] + tau * second_derivative(tau, coord[i], velocity[i]);
		coord[i + 1] = coord[i] + tau * velocity[i];
		/*
		double step = i * tau;
		std::cout << "coord[" << step << "] = " << coord[i] 
		<< ",\tvelocity[" << step << "] = " << velocity[i] << std::endl;
		*/
	}
}

void RungeKutta(double tau, double * coord, double * velocity){
	coord[0]    = U0;
	velocity[0] = V0;
	for(int i(0); i < PERIOD_NUMBER * N; ++i)
	{
		double k1 = tau * second_derivative(tau, coord[i], velocity[i]);
		double k2 = tau * second_derivative(tau, coord[i] + tau / 2.0, velocity[i] + k1 / 2.0);
		double k3 = tau * second_derivative(tau, coord[i] + tau / 2.0, velocity[i] + k2 / 2.0);
		double k4 = tau * second_derivative(tau, coord[i] + tau, velocity[i] + k3);
		velocity[i + 1] = velocity[i] + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
		coord[i + 1] = coord[i] + tau * velocity[i];
	}
	
}

void Print_data(double tau, double * coordE, double * velocityE, double * coordRK, double * velocityRK, double * coordS, double * velocityS)
{
	for(int i(0); i <= PERIOD_NUMBER * N; ++i)
	{
		std::cout << tau * i << ",\t" << coordE[i] << ",\t" << coordRK[i] << ",\t" << coordS[i] << std::endl;
		fflush(stdout);
		// fout << i * tau << ";" << coord[i] << std::endl;
	}
}

void Simple_solution(double tau, double* coordS)
{
	double OMEGA = sqrt(double(G) / double(L) - (double(LAMBDA) / double(L)) * (double(LAMBDA) / double(L)));
	for(int i(0); i <= PERIOD_NUMBER * N; ++i)
	{
		double A = double(V0) / double(OMEGA);
		double coord = A * exp(double(-LAMBDA * i * tau) / double(L)) * sin(double(OMEGA * i * tau));
		coordS[i] = coord;
		// fout << "coord[" << i << "] = " << coord << std::endl;
		//fout << i * tau << ";" << coord << std::endl;
	}
}




int main()
{


	// tau init
	std::cout << "Hello!\n";
	double Period_T = 2.0 * M_PI * sqrt((double(L) / double(G)));
	std::cout << "Period_T = " <<  Period_T << std::endl;
	
	double tau = Period_T / N;
	std::cout << "tau = " <<  tau << std::endl;

	// array mem alloc for N + 1 points from one period for N tau
	double* coordE =    new double[PERIOD_NUMBER * N + 1];
	double* velocityE = new double[PERIOD_NUMBER * N + 1];
	double* coordRK =    new double[PERIOD_NUMBER * N + 1];
	double* velocityRK = new double[PERIOD_NUMBER * N + 1];
	double* coordS =    new double[PERIOD_NUMBER * N + 1];
	double* velocityS = new double[PERIOD_NUMBER * N + 1];

	EulerSolution(tau, coordE, velocityE);
	RungeKutta(tau, coordRK, velocityRK);
	Simple_solution(tau, coordS);

	Print_data(tau, coordE, velocityE, coordRK, velocityRK, coordS, velocityS);

	delete[] coordE;
	delete[] velocityE;
	delete[] coordRK;
	delete[] velocityRK;
	delete[] coordS;
	delete[] velocityS;
}