#ifndef __MONTECARLOENGINE_H
#define __MONTECARLOENGINE_H

#include <boost/function.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/accumulators/statistics/variance.hpp>


namespace tarnpricing {

template<typename T>
struct BoostAccumulator
{
	typedef boost::accumulators::accumulator_set<T, boost::accumulators::stats<boost::accumulators::tag::variance(boost::accumulators::lazy)> > type;
};

template<typename ARG, typename RETTYPE, typename RNG, typename ACCUMTYPE = BoostAccumulator<RETTYPE>::type>
class MonteCarloEngine {
	// handy type definitions
	typedef boost::function<RETTYPE (ARG)> Functor;

	// boost::math::normal normDist;
	Functor functor;
	RNG rng; // random number generator

public:
	// constructor
	MonteCarloEngine(Functor a_functor, RNG& a_rng): functor(a_functor), rng(a_rng) {}

	// api methods
	void simulate(ACCUMTYPE& accumulator, unsigned long numSimulations);
	// void sumulate(ACCUMTYPE& accumulator, unsigned long initialNumSimulations, double requiredAccuracy, double confidence = .95);
};

template<typename ARG, typename RETTYPE, typename RNG, typename ACCUMTYPE>
void MonteCarloEngine<ARG, RETTYPE, RNG, ACCUMTYPE>::simulate(ACCUMTYPE& accumulator, unsigned long numSimulations)
{
	for (unsigned long i = 0; i < numSimulations; i++)
		accumulator(functor(rng.random()));
}

} // namespace tarnpricing

#endif