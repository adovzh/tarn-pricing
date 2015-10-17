#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <blitz/array.h>
#include <random/normal.h>
#include "InterfaceCLAPACK.hpp"

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
	
	if (!v.empty())
	{
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	else out << "[]";


	return out;
}

typedef blitz::Array<double,1> Vector;
typedef blitz::Array<double,2> Matrix;

namespace tensor {
	blitz::firstIndex i;
	blitz::secondIndex j;
	blitz::thirdIndex k;
}

void m_diag(const Vector& vec, Matrix& m)
{
	int dim = vec.extent(blitz::firstDim);
	m.resize(dim, dim);
	m = 0.0;

	for (int i = 0; i < dim; i++) m(i, i) = vec(i);
}

void m_mult(const Matrix& m1, const Matrix& m2, Matrix& out)
{
	int m1_d1 = m1.extent(blitz::firstDim);
	int m1_d2 = m1.extent(blitz::secondDim);
	int m2_d1 = m2.extent(blitz::firstDim);
	int m2_d2 = m2.extent(blitz::secondDim);

	if (m1_d2 != m2_d1)
	{
		std::stringstream error_message;
		error_message << "Matrices are inconsistent: (" << m1_d1 << ", " << m1_d2 << ") x (" << m2_d1 << "," << m2_d2 << ")";
		throw(std::logic_error(error_message.str()));
	}

	out.resize(m1_d1, m2_d2);
	out = sum(m1(tensor::i, tensor::k) * m2(tensor::k, tensor::j), tensor::k);
}

int main()
{
	const int N = 5;
	const int freq = 1;
	const int dim = N * freq;
	Matrix rho(dim, dim);


	rho = exp(-.5 * abs((tensor::i - tensor::j) / double(freq)));

	// std::cout << rho << std::endl;

	Vector eigenval;
	Matrix eigenvec;

	quantfin::interfaceCLAPACK::SymmetricEigenvalueProblem(rho, eigenval, eigenvec);

	std::cout << "Eigenval: " << eigenval << std::endl;
	std::cout << "Eigenvec: " << eigenvec << std::endl;

	Vector lambda(eigenval.extent(blitz::firstDim));
	
	lambda = sqrt(eigenval);
	std::cout << "Lambda: " << lambda << std::endl;

	Vector sorted_lambda(dim);
	blitz::Range all = blitz::Range::all();
	blitz::Range desc = blitz::Range(lambda.extent(blitz::firstDim) - 1, 0, -1);
	sorted_lambda = lambda(desc);

	std::cout << "Desc Lambda: " << sorted_lambda << std::endl;

	Matrix desc_eigenvec = eigenvec(all, desc);
	std::cout << "Desc EigenVec: " << desc_eigenvec << std::endl;

	Matrix diag_lambda;
	m_diag(lambda, diag_lambda);

	std::cout << "Diag Lambda: " << diag_lambda << std::endl;

	Matrix m1(2,5);
	m1 = 1,2,3,4,5,6,7,8,9,10;
	std::cout << "m1" << m1 << std::endl;

	Matrix m2(5,3);
	m2 = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15;
	std::cout << "m2" << m2 << std::endl;

	Matrix m3;
	m_mult(m1, m2, m3);
	std::cout << "m1 x m2 " << m3 << std::endl;
}

//int main()
//{
//	std::cout << "Hello, World!" << std::endl;
//
//	std::vector<double> v1;
//	v1.reserve(15);
//	std::cout << "Vector capacity is " << v1.capacity() << std::endl;
//	std::cout << "Vector size is " << v1.size() << std::endl;
//	std::cout << "v1: " << v1 << std::endl;
//
//	for (int i = 1; i <= 15; i++)
//	{
//		v1.push_back(double(i) / 15);
//	}
//
//	std::cout << "v1: " << v1 << std::endl;
//
//	blitz::Array<double,1> A(15);
//	blitz::firstIndex i;
//	A = (i + 1.0) / 15;
//	std::cout << "A" << A << std::endl;
//
//	blitz::Array<double,1> B(15);
//	ranlib::NormalUnit<> normal_rng;
//	for (blitz::Array<double,1>::iterator it = B.begin(); it != B.end(); it++) *it = normal_rng.random();
//	std::cout << "B" << B << std::endl;
//}