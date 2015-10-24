#include <iostream>
#include "LIBORMarketModel.h"

namespace tarnpricing {

void LIBORMarketModel::operator()(const RealMatrix& x, LowerTriangularMatrix& underlyingValues)
{
	std::cout << "Number of points: " << timeline->length() << std::endl;
}

}