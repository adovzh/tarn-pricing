#ifndef __TARNPRICING_MONTECARLOENGINE_H
#define __TARNPRICING_MONTECARLOENGINE_H

#include <boost/function.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "Types.h"


namespace tarnpricing {

/**
 * Boost accumulator of type T. Used to gather results from Monte-Carlo simulation.
 */
template<typename T>
struct BoostAccumulator
{
	typedef boost::accumulators::accumulator_set<T, boost::accumulators::stats<boost::accumulators::tag::variance(boost::accumulators::lazy)> > type;
};

/**
 * Generic Monte-Carlo engine implementation. Template parameter ARG specifies the type of the random draw per simulation,
 * RETTYPE specifies the type of the payoff, RNG specifies a functor type, used to generate a random draw (defaults to Boost 
 * functor) and ACCUMTYPE specifies the type of the accumulator (defaults to Boost accumulator).
 */
template<typename ARG, typename RETTYPE, typename RNG = BoostRNG<ARG>::type, typename ACCUMTYPE = BoostAccumulator<RETTYPE>::type>
class MonteCarloEngine {
	// handy type definitions
	typedef boost::function<RETTYPE (ARG)> Functor;

	Functor functor; ///< functor, specifies the mapping between random draw and a payoff
	RNG rng; ///< random number generator

public:
	// constructor
	MonteCarloEngine(Functor a_functor, RNG& a_rng): functor(a_functor), rng(a_rng) {}

	// simulate numSimulations times and accumulate the results
	void simulate(ACCUMTYPE& accumulator, unsigned long numSimulations);
};

template<typename ARG, typename RETTYPE, typename RNG, typename ACCUMTYPE>
void MonteCarloEngine<ARG, RETTYPE, RNG, ACCUMTYPE>::simulate(ACCUMTYPE& accumulator, unsigned long numSimulations)
{
	// rand -> functor -> accumulator
	for (unsigned long i = 0; i < numSimulations; i++)
		accumulator(functor(rng()));
}

} // namespace tarnpricing

#endif // __TARNPRICING_MONTECARLOENGINE_H