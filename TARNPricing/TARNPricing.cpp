#include <iostream>
#include <random/normal.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions/normal.hpp>

#include "Logging.h"
#include "MonteCarloEngine.h"
#include "RandomMatrix.h"
#include "Timeline.h"
#include "LowerTriangularMatrix.h"
#include "LIBORMarketModel.h"
#include "ParameterisedVolatility.h"
#include "Mapping.h"
#include "CapletPayoff.h"

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
	try 
	{
		const int dim = 3;
		const int n = 2;

		RealVector rvPoints(n + 1);
		rvPoints = 0, 10, 10.5;
		Timeline::ConstPtr timeline = Timeline::ConstPtr(new Timeline(rvPoints));

		ranlib::NormalUnit<double> normal;
		// normal.seed((unsigned int)time(0));
		BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &normal);
		RandomMatrix<double> matrix_rng(rng, dim, n);

		RealVector initial(n);
		initial = 0.05;

		RealVector a(dim), c(dim), v(dim);
		a = 0.1; c = 0.1;
		ParameterisedVolatility::ConstPtr pvol(new ParameterisedVolatility(dim, timeline, a, c));

		LIBORMarketModel<ParameterisedVolatility>::Ptr model(new LIBORMarketModel<ParameterisedVolatility>(timeline, pvol, initial));
		CapletPayoff::ConstPtr payoff(new CapletPayoff(0.0, 10.0, 10.5, .05));
		Mapping<LIBORMarketModel<ParameterisedVolatility>,RealMatrix> mapping(model, payoff);

		MonteCarloEngine<RealMatrix, double, RandomMatrix<double> > engine(boost::bind(&Mapping<LIBORMarketModel<ParameterisedVolatility>,RealMatrix>::mapping, &mapping, _1), matrix_rng);
		BoostAccumulator<double>::type accumulator;

		engine.simulate(accumulator, 10000);
	
		INFO_MESSAGE("Mean: "<< mean(accumulator))
		INFO_MESSAGE("Stdev: " << std::sqrt(variance(accumulator) / (count(accumulator) - 1)))

		boost::math::normal s;
		RealVector sigmaVec(dim);
		(*pvol)(0, 1, sigmaVec);
		double sigma2 = blitz::sum(sigmaVec * sigmaVec);
		double d1 = (log(.05 / .05) + sigma2 * 10 / 2) / sqrt(sigma2 * 10);
		double d2 = (log(.05 / .05) - sigma2 * 10 / 2) / sqrt(sigma2 * 10);
		double x = (.05 * boost::math::cdf(s, d1) - .05 * boost::math::cdf(s, d2)) * .5 / ((1 + 10 * .05) * (1 + 0.5 * 0.05));
		INFO_MESSAGE("theoretical: " << x)
	}
	catch (std::logic_error e)
	{
		std::cerr << e.what() << std::endl;
	}

	return 0;
}

//int main4()
//{
//	try
//	{
//		const int n = 1;
//		const int dim = 1;
//
//		// Timeline::Ptr timeline = Timeline::Ptr(new Timeline(0, 10, n));
//		RealVector rvPoints(3);
//		rvPoints = 0, 10, 10.5;
//		Timeline::Ptr timeline = Timeline::Ptr(new Timeline(rvPoints));
//
//		BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &ranlib::NormalUnit<double>());
//		RandomMatrix<double> matrix_rng(rng, dim, n);
//		RealMatrix randomDraw = matrix_rng();
//		LowerTriangularMatrix rates(timeline->length());
//
//		RealVector initial(timeline->length());
//		initial = 0.05;
//		rates.setColumn(0, initial);
//
//		RealVector a(dim), c(dim), v(dim);
//		a = 0.1; c = 0.1;
//		ParameterisedVolatility::ConstPtr pvol(new ParameterisedVolatility (dim, timeline, a, c));
//
//		DEBUG_MESSAGE("v = " << v)
//
//		LIBORMarketModel<ParameterisedVolatility> model(timeline, pvol);
//
//#ifdef INFO_ENABLED
//		for (int i = 0; i < timeline->length(); i++) {
//			std::cout << __LOG_PREFIX__;
//			for (int j = 0; j <= i; j++) {
//				std::cout << rates(i, j) << " ";
//			}
//			std::cout << std::endl;
//		}
//#endif
//
//		model(randomDraw, rates);
//
//#ifdef INFO_ENABLED
//		for (int i = 0; i < timeline->length(); i++) {
//			std::cout << __LOG_PREFIX__;
//			for (int j = 0; j <= i; j++) {
//				std::cout << rates(i, j) << " ";
//			}
//			std::cout << std::endl;
//		}
//#endif
//	}
//	catch (std::logic_error e) {
//		std::cerr << e.what() << std::endl;
//	}
//	return 0;
//}

//int main3()
//{
//	Timeline::Ptr timeline = Timeline::Ptr(new Timeline(0, 5, 10));
//	LIBORMarketModel model(timeline);
//
//	BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &ranlib::NormalUnit<double>());
//	RandomMatrix<double> matrix_rng(rng, 3, 5);
//	RealMatrix randomDraw = matrix_rng();
//	LowerTriangularMatrix rates(timeline->length());
//
//	int shift = 4;
//	RealVector initial(timeline->length() - shift);
//	initial = 0.05;
//	rates.setColumn(shift, initial);
//
//	std::cout << "Underlying rates" << rates << std::endl;
//
//	model(randomDraw, rates);
//
//	for (int i = 0; i < timeline->length(); i++) {
//		for (int j = 0; j <= i; j++) {
//			std::cout << rates(i, j) << " ";
//		}
//		std::cout << std::endl;
//	}
//
//	for (int i = 0; i < timeline->length(); i++)
//	{
//		RealVector v(timeline->length() - i);
//		rates.getColumn(i, v);
//		std::cout << "Column " << i << ": " << v << std::endl;
//	}
//
//	return 0;
//}
//
//int main2()
//{
//	BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &ranlib::NormalUnit<double>());
//	RandomMatrix<double> matrix_rng(rng, 3, 5);
//
//	Matrix<double>::type& matrix = matrix_rng();
//	std::cout << matrix << std::endl;
//
//	Timeline T(0.0, 5.0, 10);
//	std::cout << "Timeline: " << T << std::endl;
//
//	blitz::Array<double,1> time_points(5);
//	time_points = .35, .85, 1.15, 1.75, 1.89;
//	Timeline T2(time_points);
//	std::cout << "Timeline: " << T2 << std::endl;
//
//	LowerTriangularMatrix m(4);
//
//	for (int i = 0, s = 0; i < 4; i++)
//		for (int j = 0; j <= i; j++, s++)
//			m(i, j) = s + .07;
//
//	// std::cout << "Underlying vector: " << m.values << std::endl;
//
//	for (int i = 0; i < 4; i++)
//		for (int j = 0; j <= i; j++)
//			std::cout << "m(" << i << "," << j << ") = " << m(i, j) << std::endl;
//
//	return 0;
//}
//
//int main1()
//{
//	// ranlib::NormalUnit<double> normalDist;
//	BoostRNG<double>::type rng = boost::bind(std::mem_fun(&ranlib::NormalUnit<double>::random), &ranlib::NormalUnit<double>());
//	BoostAccumulator<double>::type accumulator;
//	MonteCarloEngine<double,double> engine(identity, rng);
//
//	const unsigned long nSimulations = 10000;
//	engine.simulate(accumulator, nSimulations);
//
//	std::cout << "Hello, World !" << std::endl;
//	std::cout << "Result: " << mean(accumulator) << std::endl;
//	
//	SimpleAccum testAccum;
//	accumulator_set<double, stats<tag::variance(lazy) > > boostAccumulator;
//
//	for (double d = 1.2, i = 0; i < 10; d += 1.1, i++)
//	{
//		testAccum(d);
//		boostAccumulator(d);
//	}
//
//	std::cout << "Mean: " << testAccum.mean() << "\t" << mean(boostAccumulator) << std::endl;
//	std::cout << "Stdev: " << testAccum.stdev() << "\t" << std::sqrt(variance(boostAccumulator) / (count(boostAccumulator) - 1)) << std::endl;
//
//	return 0;
//}