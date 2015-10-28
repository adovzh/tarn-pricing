#ifndef __LIBOR_MARKET_MODEL_H
#define __LIBOR_MARKET_MODEL_H

#include <boost/shared_ptr.hpp>
#include "Logging.h"
#include "Types.h"
#include "Timeline.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

/**
 * Template that simulates initial term structure of forward LIBOR rates through across the timeline.
 * For simplicity, the same timeline is used for both term structure of LIBOR rates and their evolution.
 *
 * The template parameter VOL class must implement a member function dimension() 
 * and define operator()(int i, int j, RealVector& sigma), that fills a vector sigma with the values of the volatility
 * at time T_i for the rate L_j. ParameterisedVolatility class satisfied these requirements.
 */
template<typename VOL, typename VOLCPTR = VOL::ConstPtr>
class LIBORMarketModel
{
public:
	LIBORMarketModel(const Timeline::ConstPtr& _timeline, const VOLCPTR& _vol, const RealVector& _initial): m_timeline(_timeline), pvol(_vol), m_initial(_initial) {}
	/// Generate the realisation
	void operator()(const RealMatrix& x, LowerTriangularMatrix& underlyingValues);
	/// retrieves the timeline
	const Timeline::ConstPtr& timeline() const { return m_timeline; }
	/// sets the timeline
	void timeline(const Timeline::ConstPtr& _timeline) { m_timeline = _timeline; }
	/// retrieves the initial
	const RealVector& initial() const { return m_initial; }
	/// sets initial vector
	void initial(const RealVector& _initial) { m_initial = _initial; }

	/// type definitions
	typedef boost::shared_ptr<LIBORMarketModel> Ptr;
	typedef boost::shared_ptr<const LIBORMarketModel> ConstPtr;
private:
	Timeline::ConstPtr m_timeline; ///< Timeline for both term struture of simulated rates and their evolution.
	const VOLCPTR pvol;	///< Smart pointer to a constant volatility.
	RealVector m_initial; ///< Initial term structure of forward LIBOR rates.

	/// calculate drift term
	double calculateDrift(const LowerTriangularMatrix& underlying, const RealVector& sj, double vNormSq, int i, int j);
};

template <typename VOL, typename VOLCPTR>
void LIBORMarketModel<VOL, VOLCPTR>::operator()(const RealMatrix& x, LowerTriangularMatrix& underlyingValues)
{
	TRACE_MESSAGE("Number of points: " << timeline()->length())
	int n = timeline()->length();
	const VOL& vol = *pvol;

	for (int i = 0; i < n - 1; i++)
	{
		TRACE_MESSAGE("Simulation step T[" << i << "]")
		RealVector z = x(blitz::Range::all(), i);

		for (int j = i + 1; j < n; j++)
		{
			TRACE_MESSAGE("Simulating L_" << j << "(t_" << i + 1 << ')')
			RealVector v(vol.dimension());
			vol(i, j, v);

			// sigma and Z dot-product
			double sz = blitz::sum(v * z);
			TRACE_MESSAGE("sz: " << sz)

			// calculate drift
			// pre-calculate sigma norm squared
			double vNormSq = blitz::sum(v * v);
			double drift = calculateDrift(underlyingValues, v, vNormSq, i, j);

			double delta = timeline()->delta(i);
			TRACE_MESSAGE("delta_" << i << ": " << delta)
			underlyingValues(j, i + 1) = underlyingValues(j, i) * exp((drift - vNormSq / 2) * delta + sqrt(delta) * sz);
			TRACE_MESSAGE("=========================")
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
		double deltaL = timeline()->delta(k) * underlying(k, i);
		double driftTerm = blitz::sum(sj * sk) * deltaL / (1 + deltaL);

		drift += driftTerm;
	}

	// now the last term
	deltaL = timeline()->delta(j) * underlying(j, i);
	driftTerm = vNormSq * deltaL / (1 + deltaL);

	drift += driftTerm;

	return drift;
}

} // tarnprincing

#endif // __LIBOR_MARKET_MODEL_H