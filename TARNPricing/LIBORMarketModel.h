#ifndef __LIBOR_MARKET_MODEL_H
#define __LIBOR_MARKET_MODEL_H

#include "Types.h"
#include "Timeline.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

class LIBORMarketModel
{
private:
	const Timeline::Ptr timeline;
public:
	LIBORMarketModel(const Timeline::Ptr& _timeline): timeline(_timeline) {}
	// Generate the realisation
	void operator()(const RealMatrix& x, LowerTriangularMatrix& underlyingValues);
};

} // tarnprincing

#endif // __LIBOR_MARKET_MODEL_H