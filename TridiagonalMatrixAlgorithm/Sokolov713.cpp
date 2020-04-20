

#include <cstddef>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <climits>

static const double A = 1.0;
static const double B = 2.0 * exp(1);

class Matrix
{
public:
	Matrix(){};
	Matrix(std::size_t dimension, const std::string & filename);
	~Matrix();
	std::size_t GetDim() const;
	void PrintVector();
	void SaveResults(std::vector<double> & approximate_solution);
	std::vector<double> TridiagonalAlgo();
	double CalculateNorm(std::vector<double> & approximate_solution);
private:
	std::size_t dimension_;
	std::vector<double> up_;
	std::vector<double> middle_;
	std::vector<double> down_;
	std::ofstream fout_;
};

bool TridiagonalAlgo(Matrix & matrix);
long int ReadArg(char * str);

int main(int argc, char** argv)
{
	if(argc < 2){
		std::cout << "Error! Number of grid intervals required!\n";
		exit(EXIT_FAILURE);
	}
	std::size_t n = ReadArg(argv[1]);
	if(n < 2){
		std::cout << "Error! Number of grid intervals should not be less than 2\n";
		exit(EXIT_FAILURE);
	}
	
	Matrix matrix(n, "output.txt");
	std::vector<double> approximate_solution = matrix.TridiagonalAlgo();
	matrix.SaveResults(approximate_solution);
	double norm = matrix.CalculateNorm(approximate_solution);
	std::cout << "Norm = " << norm << std::endl;

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


Matrix::Matrix(std::size_t dimension, const std::string & filename):
	dimension_(dimension),
	up_(dimension),
	middle_(dimension + 1),
	down_(dimension),
	fout_(filename)
{
	double h = 1.0 / (dimension_);

	for(std::size_t i(1); i < dimension_; ++i)
	{
		down_.at(i - 1) = -1.0 / (h * h);
		middle_.at(i) = 1.0 / (h * h) + 1.0 / (h * h) + 1.0;
		up_.at(i) = -1.0 / (h * h);
	}

	up_.at(0) = -2.0 / (h * h);
	middle_.at(0) = 2.0 / (h * h) + 1.0;

	middle_.at(dimension_) = 2.0 / (h * h) + 1.0;
	down_.at(dimension_ - 1) = -2.0 / (h * h);
	
}


Matrix::~Matrix()
{
	up_.clear();
	middle_.clear();
	down_.clear();
	fout_.close();
}


void Matrix::PrintVector()
{
	std::cout << "==============================================\n";
	for (auto i : up_)
	   std::cout << i << ' ';
	std::cout << std::endl;

	for (auto i : middle_)
	   std::cout << i << ' ';
	std::cout << std::endl;

	for (auto i : down_)
	   std::cout << i << ' ';
	std::cout << std::endl;
	std::cout << "==============================================\n";
}

std::size_t Matrix::GetDim() const
{
	return this->dimension_;
}

double function(double x)
{
	return exp(x);
}

double exact_solution(double x)
{
	return x * exp(x);
}

std::vector<double> CalculateRightSide(std::size_t dim)
{
	std::vector<double> right;
	double h = 1.0 / (dim);
	right.push_back(-2.0 - 2.0 * A / h);
	for(std::size_t i(1); i < dim; ++i){
		right.push_back(-2.0 * function(double(i) * h));
	}
	right.push_back(-2.0 * function(1) + 2.0 * B / h);

	return right;
}

double Matrix::CalculateNorm(std::vector<double> & approximate_solution)
{
	double norm(0.0);
	for(std::size_t i(0); i < dimension_ + 1; ++i)
	{
		double x(double(i) / dimension_);
		double delta = fabs(approximate_solution[i] - exact_solution(x));
		if(delta > norm)
			norm = delta;
	}
	return norm;
}


std::vector<double> Matrix::TridiagonalAlgo()
{
	std::vector<double> right = CalculateRightSide(dimension_);
    
    up_.at(0) /= middle_.at(0);
    right.at(0) /= middle_.at(0);

    for(std::size_t i(1); i < dimension_; ++i){
    	up_.at(i) /= middle_.at(i) - down_.at(i-1) * up_.at(i - 1);
    	right.at(i) = (right.at(i) - down_.at(i-1) * right.at(i - 1)) / (middle_.at(i) - down_.at(i-1) * up_.at(i - 1));
    }
    right.at(dimension_) = (right.at(dimension_) - down_.at(dimension_ - 1) * right.at(dimension_ - 1)) / (middle_.at(dimension_) - down_.at(dimension_ - 1) * up_.at(dimension_ - 1));
    

    for(std::size_t i = dimension_; i-- >0;){
    	right.at(i) -= up_.at(i) * right.at(i + 1);
    }

    return right;
}

void Matrix::SaveResults(std::vector<double> & approximate_solution)
{
	for(std::size_t i(0); i < dimension_ + 1; ++i)
	{
		double x(double(i) / dimension_);
		fout_ << x << " " << approximate_solution[i] << " " << exact_solution(x) << std::endl;
	}
}


