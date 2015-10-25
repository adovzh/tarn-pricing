#ifndef __LIBOR_MARKET_MODEL_H
#define __LIBOR_MARKET_MODEL_H

#include "Logging.h"
#include "Types.h"
#include "Timeline.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

template<typename VOL, typename VOLCPTR = VOL::ConstPtr>
class LIBORMarketModel
{
public:
	LIBORMarketModel(const Timeline::Ptr& _timeline, const VOLCPTR& _vol): timeline(_timeline), pvol(_vol) {}
	// Generate the realisation
	void operator()(const RealMatrix& x, LowerTriangularMatrix& underlyingValues);
private:
	const Timeline::Ptr timeline;
	const VOLCPTR pvol;

	// calculate drift term
	double calculateDrift(const LowerTriangularMatrix& underlying, const RealVector& sj, double vNormSq, int i, int j);
};

template <typename VOL, typename VOLCPTR>
void LIBORMarketModel<VOL, VOLCPTR>::operator()(const RealMatrix& x, LowerTriangularMatrix& underlyingValues)
{
	LOG_MESSAGE("Number of points: " << timeline->length())
	int n = timeline->length();
	const VOL& vol = *pvol;

	for (int i = 0; i < n - 1; i++)
	{
		LOG_MESSAGE("Simulation step T[" << i << "]")
		RealVector z = x(blitz::Range::all(), i);

		for (int j = i + 1; j < n; j++)
		{
			LOG_MESSAGE("Simulating L_" << j << "(t_" << i + 1 << ')')
			RealVector v(vol.dimension());
			vol(i, j, v);

			// sigma and Z dot-product
			double sz = blitz::sum(v * z);
			LOG_MESSAGE("sz: " << sz)

			// calculate drift
			// pre-calculate sigma norm squared
			double vNormSq = blitz::sum(v * v);
			double drift = calculateDrift(underlyingValues, v, vNormSq, i, j);

			double delta = timeline->delta(i);
			underlyingValues(j, i + 1) = underlyingValues(j, i) * exp((drift - vNormSq / 2) * delta + sqrt(delta) * sz);
			LOG_MESSAGE("=========================")
		}
	}
}

template <typename VOL, typename VOLCPTR>
inline double LIBORMarketModel<VOL, VOLCPTR>::calculateDrift(const LowerTriangularMatrix& underlying, 
													const RealVector& sj, double vNormSq, 
													int i, int j)
{
	const VOL& vol = *pvol;
	double drift = 0.0;
	double deltaL, driftTerm;

	// summing up terms of a drift, except for the last one, k = j
	for (int k = i + 1; k < j; k++)
	{
		RealVector sk(vol.dimension());
		vol(i, k, sk);

		// pre-calculating delta * L
		double deltaL = timeline->delta(k) * underlying(k, i);
		double driftTerm = blitz::sum(sj * sk) * deltaL / (1 + deltaL);

		drift += driftTerm;
	}

	// now the last term
	deltaL = timeline->delta(j) * underlying(j, i);
	// driftTerm = blitz::sum(sj * sj) * deltaL / (1 + deltaL);
	driftTerm = vNormSq * deltaL / (1 + deltaL);

	drift += driftTerm;

	return drift;
}

} // tarnprincing

#endif // __LIBOR_MARKET_MODEL_H