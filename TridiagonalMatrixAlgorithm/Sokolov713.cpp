#include <cstddef>
#include <vector>
#include <iostream>
#include <climits>
#include <cmath>

class Matrix
{
public:
	Matrix(){};
	Matrix(std::size_t dimension);
	~Matrix();
	std::size_t GetDim() const;
	void PrintVector();
	bool TridiagonalAlgo();
private:
	std::size_t dimension_;
	std::vector<double> up_;
	std::vector<double> middle_;
	std::vector<double> down_;
};

bool TridiagonalAlgo(Matrix & matrix);


Matrix::Matrix(std::size_t dimension):
	dimension_(dimension)
{
	double h = 1.0 / (dimension_ - 1);

	// fill up

	// up_.push_back(0);
	for(std::size_t i(0); i < dimension_ - 1; ++i)
		up_.push_back(1 - h / 2.0);
	up_.push_back(0);

	up_.at(0) = 2.0 * h;

	// fill middle

	// middle_.push_back(1);
	for(std::size_t i(0); i < dimension_ ; ++i)
		middle_.push_back(-2);

	middle_.at(0) =  -2.0 * h;
	middle_.at(dimension_ - 1) = h;

	// fill down
	down_.push_back(0);
	for(std::size_t i(0); i < dimension_ - 1; ++i)
		down_.push_back(1 + h / 2.0);

	down_.at(dimension_ - 1) = -h;
	
}

void Matrix::PrintVector()
{
	for (auto i : up_)
	   std::cout << i << ' ';
	std::cout << std::endl;

	for (auto i : middle_)
	   std::cout << i << ' ';
	std::cout << std::endl;

	for (auto i : down_)
	   std::cout << i << ' ';
	std::cout << std::endl;
}

Matrix::~Matrix()
{
	up_.clear();
	middle_.clear();
	down_.clear();
}

std::size_t Matrix::GetDim() const
{
	return this->dimension_;
}

double function(double x)
{
	return exp(x);
}

std::vector<double> CalculateRightSide(std::size_t dim)
{
	std::vector<double> right;
	for(std::size_t i(0); i < dim; ++i){
		double h = 1.0 / (dim - 1);
		right.push_back(2.0 * h * h * exp(double(i) * h));
	}

	return right;
}


bool Matrix::TridiagonalAlgo()
{
	this->PrintVector();
	
	
	std::vector<double> right = CalculateRightSide(dimension_);
	
	/*
	std::cout << "Right side:\n";
	for (auto i : right)
	   std::cout << i << ' ';
	std::cout << std::endl;
	*/

	/*std::vector<double> new_up;
	std::vector<double> new_right;
	std::vector<double> y;

	new_up.push_back(up_.front() / middle_.front());
	new_right.push_back(right.front() / middle_.front());
	
	for(std::size_t i(1); i < dimension_ - 1; ++i)
	{
		new_up.push_back(up_.at(i) / (middle_.at(i) - down_.at(i) * new_up.at(i - 1)));
		new_right.push_back((right.at(i) - down_.at(i) * new_right.at(i - 1)) / (middle_.at(i) - down_.at(i) * new_up.at(i - 1)));
	}

	new_right.push_back((right.at(dimension_ - 1) - down_.at(dimension_ - 1) * new_right.at(dimension_ - 2)) / (middle_.at(dimension_ - 1) - down_.at(dimension_ - 1) * new_up.at(dimension_ - 2)));

	std::cout << "new_up:\n";
	for (auto i : new_up)
	   std::cout << i << '\t';
	std::cout << std::endl;

	std::cout << "new_right:\n";
	for (auto i : new_right)
	   std::cout << i << '\t';
	std::cout << std::endl;


	y.resize(dimension_);
	y.at(dimension_ - 1) = new_right.at(dimension_ - 1);

	for(int i(dimension_ - 2); i >= 0; --i)
	{
		y.at(i) = new_right.at(i) - new_up.at(i) * y.at(i + 1);
	}

	std::cout << "y\n";
	for (auto i : y)
	   std::cout << i << ' ';
	std::cout << std::endl;*/


	/*n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }*/

    dimension_--;
    
    up_.at(0) /= middle_.at(0);
    right.at(0) /= middle_.at(0);

    for(std::size_t i(1); i < dimension_; ++i){
    	up_.at(i) /= middle_.at(i) - down_.at(i) * up_.at(i - 1);
    	right.at(i) = (right.at(i) - down_.at(i) * right.at(i - 1)) / (middle_.at(i) - down_.at(i) * up_.at(i - 1));
    }
    right.at(dimension_) = (right.at(dimension_) - down_.at(dimension_) * right.at(dimension_ - 1)) / (middle_.at(dimension_) - down_.at(dimension_) * up_.at(dimension_ - 1));
    

    for(std::size_t i = dimension_; i-- >0;){
    	right.at(i) -= up_.at(i) * right.at(i + 1);
    }

    std::cout << "y\n";
	for (auto i : right)
	   std::cout << i << ' ';
	std::cout << std::endl;

	return true;
}


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
	
	Matrix matrix(n);
	matrix.TridiagonalAlgo();

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