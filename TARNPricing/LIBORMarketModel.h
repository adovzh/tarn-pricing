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
private:
	const Timeline::Ptr timeline;
	VOLCPTR pvol;
public:
	LIBORMarketModel(const Timeline::Ptr& _timeline, const VOLCPTR& _vol): timeline(_timeline), pvol(_vol) {}
	// Generate the realisation
	void operator()(const RealMatrix& x, LowerTriangularMatrix& underlyingValues);
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

			double sz = blitz::sum(v * z);
			LOG_MESSAGE("sz: " << sz)

			double delta = timeline->delta(i);
			underlyingValues(j, i + 1) = underlyingValues(j, i) * exp(sqrt(delta) * sz);
			LOG_MESSAGE("=========================")
		}
	}
}


} // tarnprincing

#endif // __LIBOR_MARKET_MODEL_H