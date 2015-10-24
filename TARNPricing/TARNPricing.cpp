#include <iostream>
#include <random/normal.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include "MonteCarloEngine.h"
#include "RandomMatrix.h"
#include "Timeline.h"
#include "LowerTriangularMatrix.h"
#include "LIBORMarketModel.h"

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
	Timeline::Ptr timeline = Timeline::Ptr(new Timeline(0, 5, 10));
	LIBORMarketModel model(timeline);

	BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &ranlib::NormalUnit<double>());
	RandomMatrix<double> matrix_rng(rng, 3, 5);
	RealMatrix randomDraw = matrix_rng();
	LowerTriangularMatrix rates(timeline->length());

	int shift = 10;
	RealVector initial(timeline->length() - shift);
	initial = 0.05;
	rates.setColumn(shift, initial);

	std::cout << "Underlying rates" << rates << std::endl;

	model(randomDraw, rates);

	for (int i = 0; i < timeline->length(); i++) {
		for (int j = 0; j <= i; j++) {
			std::cout << rates(i, j) << " ";
		}
		std::cout << std::endl;
	}

	return 0;
}

int main2()
{
	BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &ranlib::NormalUnit<double>());
	RandomMatrix<double> matrix_rng(rng, 3, 5);

	Matrix<double>::type& matrix = matrix_rng();
	std::cout << matrix << std::endl;

	Timeline T(0.0, 5.0, 10);
	std::cout << "Timeline: " << T << std::endl;

	blitz::Array<double,1> time_points(5);
	time_points = .35, .85, 1.15, 1.75, 1.89;
	Timeline T2(time_points);
	std::cout << "Timeline: " << T2 << std::endl;

	LowerTriangularMatrix m(4);

	for (int i = 0, s = 0; i < 4; i++)
		for (int j = 0; j <= i; j++, s++)
			m(i, j) = s + .07;

	// std::cout << "Underlying vector: " << m.values << std::endl;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j <= i; j++)
			std::cout << "m(" << i << "," << j << ") = " << m(i, j) << std::endl;

	return 0;
}

int main1()
{
	// ranlib::NormalUnit<double> normalDist;
	BoostRNG<double>::type rng = boost::bind(std::mem_fun(&ranlib::NormalUnit<double>::random), &ranlib::NormalUnit<double>());
	BoostAccumulator<double>::type accumulator;
	MonteCarloEngine<double,double> engine(identity, rng);

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

	return 0;
}