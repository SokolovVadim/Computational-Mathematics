#include <cstddef>
#include <vector>

class Matrix
{
public:
	Matrix(){};
	Matrix(std::size_t dimension);
	~Matrix();
	std::size_t GetDim() const;
	void PrintVector();
private:
	std::size_t dimension_;
	std::vector<double> up_;
	std::vector<double> middle_;
	std::vector<double> down_;
};

bool TridiagonalAlgo(Matrix & matrix);