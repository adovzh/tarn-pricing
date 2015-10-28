#ifndef __TARNPRICING_RANDOMMATRIX_H
#define __TARNPRICING_RANDOMMATRIX_H

#include "Types.h"

namespace tarnpricing {

/**
 * Random matrix generator, inspired by RandomArray implementation.
 */
template<typename rv_type, typename RNG = BoostRNG<rv_type>::type>
class RandomMatrix
{
private:
	typename Matrix<rv_type>::type random_values;
	RNG rng;
public:
	RandomMatrix(RNG& _rng, int rows, int cols): random_values(rows, cols), rng(_rng) {}
	typename Matrix<rv_type>::type& operator()();
};

template<typename rv_type, typename RNG>
typename Matrix<rv_type>::type& RandomMatrix<rv_type, RNG>::operator()()
{
	for (int i = 0; i < random_values.extent(blitz::firstDim); i++)
		for (int j = 0; j < random_values.extent(blitz::secondDim); j++)
			random_values(i, j) = rng();

	return random_values;
}

} // namespace tarnpricing

#endif // __TARNPRICING_RANDOMMATRIX_H