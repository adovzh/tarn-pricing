#ifndef __LOWER_TRIANGULAR_MATRIX_H
#define __LOWER_TRIANGULAR_MATRIX_H

#include <iostream>
#include "Types.h"

namespace tarnpricing {

class LowerTriangularMatrix
{
private:
	int rows;
	RealVector values;
public:
	LowerTriangularMatrix(int r): rows(r), values(r * (r + 1) / 2) { values = 0.0; }
	int size() const { return rows; }
	
	double operator()(int i, int j) const;
	double& operator()(int i, int j);

	void getColumn(int j, RealVector& column);
	void setColumn(int j, const RealVector& column);

	friend std::ostream& operator<<(std::ostream& out, const LowerTriangularMatrix& timeline);
};

inline double& LowerTriangularMatrix::operator()(int i, int j)
{
	if (j > i) throw std::logic_error("j <= i condition violated");
	return values(i * (i + 1) / 2 + j);
}

inline double LowerTriangularMatrix::operator()(int i, int j) const
{
	if (j > i) throw std::logic_error("j <= i condition violated");
	return values(i * (i + 1) / 2 + j);
}

inline void LowerTriangularMatrix::getColumn(int j, RealVector& column)
{
	int i = 0, index = 0, ci = 0;

	while (i < rows)
	{
		if (i < j) { index += ++i + 1; }
		else { column(ci++) = values(index); index += ++i;}
	}		
}

inline void LowerTriangularMatrix::setColumn(int j, const RealVector& column)
{
	int i = 0, index = 0, ci = 0;

	while (i < rows)
	{
		if (i < j) { index += ++i + 1; }
		else { values(index) = column(ci++); index += ++i;}
	}		
}

} // namespace tarnpricing

#endif