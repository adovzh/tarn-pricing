#include <iostream>
#include <random/normal.h>
#include "MonteCarloEngine.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace tarnpricing;
using namespace boost::accumulators;

double dummy_functor(double arg)
{
	return 1.0;
}

double identity(double arg)
{
	return arg;
}

class SimpleAccum
{
	double running_sum;
	double sum2;
	int n;
public:
	SimpleAccum(): running_sum(0.0), sum2(0.0), n(0) {}
	SimpleAccum& operator()(double j) { running_sum += j; sum2 += j * j; n++; return *this; }
	double sum() const { return running_sum; }
	double mean() const { return running_sum / n; }
	double stdev() const { double m = mean(); return sqrt(((sum2 / n) - m * m) / (n - 1)); }
	double count() const { return n; }
};

int main()
{
	ranlib::NormalUnit<double> normalDist;
	BoostAccumulator<double>::type accumulator;
	MonteCarloEngine<double,double,ranlib::NormalUnit<double> > engine(identity, normalDist);

	const unsigned long nSimulations = 10000;
	engine.simulate(accumulator, nSimulations);

	std::cout << "Hello, World !" << std::endl;
	std::cout << "Result: " << mean(accumulator) << std::endl;
	
	SimpleAccum testAccum;
	accumulator_set<double, stats<tag::variance(lazy) > > boostAccumulator;

	for (double d = 1.2, i = 0; i < 10; d += 1.1, i++)
	{
		testAccum(d);
		boostAccumulator(d);
	}

	std::cout << "Mean: " << testAccum.mean() << "\t" << mean(boostAccumulator) << std::endl;
	std::cout << "Stdev: " << testAccum.stdev() << "\t" << std::sqrt(variance(boostAccumulator) / (count(boostAccumulator) - 1)) << std::endl;
}