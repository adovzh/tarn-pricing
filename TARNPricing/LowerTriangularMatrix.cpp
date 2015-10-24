#include "LowerTriangularMatrix.h"

namespace tarnpricing {

std::ostream& operator<<(std::ostream& out, const LowerTriangularMatrix& matrix)
{
	out << matrix.values;
	return out;
}

} // namespace tarnpricing