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
	
	double& operator()(int i, int j)
	{
		if (j > i) throw std::logic_error("j <= i condition violated");
		return values(i * (i + 1) / 2 + j);
	}

	double operator()(int i, int j) const
	{
		if (j > i) throw std::logic_error("j <= i condition violated");
		return values(i * (i + 1) / 2 + j);
	}

	void setColumn(int j, const RealVector& column)
	{
		int i = 0, index = 0, ci = 0;

		while (i < rows)
		{
			if (i < j) { index += ++i + 1; }
			else { values(index) = column(ci++); index += ++i;}
		}		
	}

	friend std::ostream& operator<<(std::ostream& out, const LowerTriangularMatrix& timeline);
};

} // namespace tarnpricing

#endif