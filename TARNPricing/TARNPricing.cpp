#include <iostream>
#include <random/normal.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#define LOG_LEVEL_INFO

#include "Logging.h"
#include "MonteCarloEngine.h"
#include "RandomMatrix.h"
#include "Timeline.h"
#include "LowerTriangularMatrix.h"
#include "LIBORMarketModel.h"
#include "ParameterisedVolatility.h"
#include "Mapping.h"
#include "CapletPayoff.h"
#include "TARNInverseFloaterPayoff.h"
#include "Config.h"

using namespace tarnpricing;

namespace {
	
	/**
	 *	Auxiliary function that copies the contents of the vector
	 *  to a newly allocated C-style array of size len.
	 *  Only first len elements of the vector will be copied.
	 *  If the size of the vector is less than len, the last element
	 *  will be replicated, so that the array always contains len elements.
	 * 
	 *	This function allocates memory, that has to be reclaimed by the caller.
	 *  One useful scenario for it is to create Blitz++ Array from the existing data.
	 *  Blitz++ will then take care of the memory.
	 */
	double* copyVectorData(int len, const std::vector<double>& v)
	{
		double* data = new double[len];
		std::vector<double>::const_iterator it = v.cbegin();
		double lastElement = 0.0;

		for (int i = 0;  i < len; i++)
			data[i] = (it != v.cend()) ? (lastElement = *it++) : lastElement;

		return data;
	}

}

namespace tarnpricing {
	// forward declarations of supported programs
	void tarn(const Config::ConstPtr&);
	void caplet(const Config::ConstPtr&);
}

int main(int argc, char* argv[])
{
	try 
	{
		// read the configuration file parameter from the command line
		char* configFile = (argc > 1) ? argv[1] : "parameters.csv";
		INFO_MESSAGE("Reading configuration from " << configFile)

		// parse configuration file
		Config::ConstPtr config = Config::parse(configFile);

		if (!config) 
		{
			std::cout << std::endl;
			std::cout << "No csv file specified or default 'parameters.csv' cannot be found." << std::endl;
			std::cout << "USAGE: TARNPricing.exe <parameters csv file>" << std::endl;
			return -1;
		}

		// read the name of the program from the config
		std::string program = config->getString("program");

		INFO_MESSAGE("Running program " << program)

		try {
			std::string description = config->getString("description");
			INFO_MESSAGE("Description: " << description)
		} catch (const std::logic_error&) {
			// ignore
		}

		if (program == "caplet")
			caplet(config);
		else if (program == "tarn")
			tarn(config);
		else
			INFO_MESSAGE("No program to run, please define 'program' in the configuration file");
	}
	catch (const boost::bad_lexical_cast& e)
	{
		ERROR_MESSAGE("Conversion error: " << e.what())
	}
	catch (std::logic_error e)
	{
		ERROR_MESSAGE(e.what())
	}

	return 0;
}

namespace tarnpricing {

// why not override opertor<<() for standard verctors?
template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
	
	if (!v.empty())
	{
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	else out << "[]";


	return out;
}

// pricing a Target Redemption Note
void tarn(const Config::ConstPtr& config)
{
	// read program-specific configuration parameters
	const int dim = config->getInt("dim");
	DEBUG_MESSAGE("Dimension of underlying Brownian Motions (dim): " << dim)

	const int n = config->getInt("n");
	DEBUG_MESSAGE("Number of LIBOR rates (n): " << n)

	std::vector<double> timeline_vec;
	config->getDoubles("timeline", timeline_vec);
	DEBUG_MESSAGE("Timeline: " << timeline_vec)

	const double strike = config->getDouble("Target");
	DEBUG_MESSAGE("Target coupon: " << strike)

	const double c1 = config->getDouble("c1");
	DEBUG_MESSAGE("c1: " << c1)

	const double K = config->getDouble("K");
	DEBUG_MESSAGE("K: " << K);

	const double alpha = config->getDouble("alpha");
	DEBUG_MESSAGE("alpha: " << alpha)

	const int nSim = config->getInt("nSim");
	DEBUG_MESSAGE("Number of simulations (nSim): " << nSim)

	std::vector<double> covar_a;
	config->getDoubles("covar_a", covar_a);
	DEBUG_MESSAGE("covar_a: " << covar_a)

	std::vector<double> covar_c;
	config->getDoubles("covar_c", covar_c);
	DEBUG_MESSAGE("covar_c: " << covar_c)

	std::vector<double> initial_rates;
	config->getDoubles("initial", initial_rates);
	DEBUG_MESSAGE("Initial forward curve (initial_rates)" << initial_rates)

	// timeline
	double* pTimeline = copyVectorData(n + 1, timeline_vec);
	RealVector rvPoints(pTimeline, blitz::shape(n + 1), blitz::deleteDataWhenDone);
	Timeline::ConstPtr timeline(new Timeline(rvPoints));
	DEBUG_MESSAGE("Time points: " << *timeline)

	// random number generator
	ranlib::NormalUnit<double> normal;
	normal.seed((unsigned int)time(0));
	BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &normal);
	RandomMatrix<double> matrix_rng(rng, dim, n);

	// converts initial LIBOR rates to Blitz++ Array
	double* initialRatesData = copyVectorData(n + 1, initial_rates);
	RealVector initial(initialRatesData, blitz::shape(n), blitz::deleteDataWhenDone);
	DEBUG_MESSAGE("initial = " << initial)

	// convert convariance parameter vectors to Blitz Arrays
	double* covar_a_data = copyVectorData(dim, covar_a);
	RealVector a(covar_a_data, blitz::shape(dim), blitz::deleteDataWhenDone);
	DEBUG_MESSAGE("a = " << a)

	double* covar_c_data = copyVectorData(dim, covar_c);
	RealVector c(covar_c_data, blitz::shape(dim), blitz::deleteDataWhenDone);
	DEBUG_MESSAGE("c = " << c)

	// volatility structure
	ParameterisedVolatility::ConstPtr pvol(new ParameterisedVolatility(dim, timeline, a, c));

	// initialise model, payoff and link them through mapping
	LIBORMarketModel<ParameterisedVolatility>::Ptr model(new LIBORMarketModel<ParameterisedVolatility>(timeline, pvol, initial));
	TARNInverseFloaterPayoff::ConstPtr payoff(new TARNInverseFloaterPayoff(timeline, strike, c1, K, alpha));
	Mapping<LIBORMarketModel<ParameterisedVolatility>,RealMatrix> mapping(model, payoff);

	// initialise Monte-Carlo simulation engine
	MonteCarloEngine<RealMatrix, double, RandomMatrix<double> > engine(boost::bind(&Mapping<LIBORMarketModel<ParameterisedVolatility>,RealMatrix>::mapping, &mapping, _1), matrix_rng);
	// results will be saved into Boost accumulator
	BoostAccumulator<double>::type accumulator;

	// run simulations
	INFO_MESSAGE("Running " << nSim << " simulations")
	engine.simulate(accumulator, nSim);
	
	// calculate statistics
	double avg = boost::accumulators::mean(accumulator);
	double stdev = std::sqrt(boost::accumulators::variance(accumulator) / (boost::accumulators::count(accumulator) - 1));
	double lowerBound = avg - 1.96 * stdev;
	double upperBound = avg + 1.96 * stdev;

	INFO_MESSAGE("Monte-Carlo estimate: " << avg)
	INFO_MESSAGE("Stdev: " << stdev)
	INFO_MESSAGE("95% confidence interval: [" << lowerBound << ", " << upperBound << "]")
} // tarnpricing::tarn

// Testing a LIBOR Market Model on a caplet
void caplet(const Config::ConstPtr& config)
{
	// read program-specific configuration parameters
	const int dim = config->getInt("dim");
	DEBUG_MESSAGE("Dimension of underlying Brownian Motions (dim): " << dim)

	const int n = config->getInt("n");
	DEBUG_MESSAGE("Number of LIBOR rates (n): " << n)

	std::vector<double> timeline_vec;
	config->getDoubles("timeline", timeline_vec);
	DEBUG_MESSAGE("Timeline: " << timeline_vec)

	const double K = config->getDouble("Strike");
	DEBUG_MESSAGE("Strike: " << K)

	const int nSim = config->getInt("nSim");
	DEBUG_MESSAGE("Number of simulations (nSim): " << nSim)

	std::vector<double> covar_a;
	config->getDoubles("covar_a", covar_a);
	DEBUG_MESSAGE("covar_a: " << covar_a)

	std::vector<double> covar_c;
	config->getDoubles("covar_c", covar_c);
	DEBUG_MESSAGE("covar_c: " << covar_c)

	std::vector<double> initial_rates;
	config->getDoubles("initial", initial_rates);
	DEBUG_MESSAGE("Initial forward curve (initial_rates)" << initial_rates)

	// timeline
	double *pTimeline = copyVectorData(n + 1, timeline_vec);
	RealVector rvPoints(pTimeline, blitz::shape(n + 1), blitz::deleteDataWhenDone);
	Timeline::ConstPtr timeline(new Timeline(rvPoints));

	// random number generator
	ranlib::NormalUnit<double> normal;
	normal.seed((unsigned int)time(0));
	BoostRNG<double>::type rng = boost::bind(&ranlib::NormalUnit<double>::random, &normal);
	RandomMatrix<double> matrix_rng(rng, dim, n);

	// converts initial LIBOR rates to Blitz++ Array
	double* initialRatesData = new double[initial_rates.size()];
	std::copy(initial_rates.begin(), initial_rates.end(), initialRatesData);
	RealVector initial(initialRatesData, blitz::shape(n), blitz::deleteDataWhenDone);
	DEBUG_MESSAGE("initial = " << initial)

	// convert convariance parameter vectors to Blitz Arrays
	double* covar_a_data = new double[covar_a.size()];
	std::copy(covar_a.begin(), covar_a.end(), covar_a_data);
	RealVector a(covar_a_data, blitz::shape(dim), blitz::deleteDataWhenDone);
	DEBUG_MESSAGE("a = " << a)

	double* covar_c_data = new double[covar_c.size()];
	std::copy(covar_c.begin(), covar_c.end(), covar_c_data);
	RealVector c(covar_c_data, blitz::shape(dim), blitz::deleteDataWhenDone);
	DEBUG_MESSAGE("c = " << c)

	// volatility structure
	ParameterisedVolatility::ConstPtr pvol(new ParameterisedVolatility(dim, timeline, a, c));

	// initialise model, payoff and link them through mapping
	LIBORMarketModel<ParameterisedVolatility>::Ptr model(new LIBORMarketModel<ParameterisedVolatility>(timeline, pvol, initial));
	CapletPayoff::ConstPtr payoff(new CapletPayoff(timeline_vec[0], timeline_vec[1], timeline_vec[2], K));
	Mapping<LIBORMarketModel<ParameterisedVolatility>,RealMatrix> mapping(model, payoff);

	// initialise Monte-Carlo simulation engine
	MonteCarloEngine<RealMatrix, double, RandomMatrix<double> > engine(boost::bind(&Mapping<LIBORMarketModel<ParameterisedVolatility>,RealMatrix>::mapping, &mapping, _1), matrix_rng);
	// results will be saved into Boost accumulator
	BoostAccumulator<double>::type accumulator;

	// run simulations
	INFO_MESSAGE("Running " << nSim << " simulations")
	engine.simulate(accumulator, nSim);
	
	// calculate statistics
	double avg = boost::accumulators::mean(accumulator);
	double stdev = std::sqrt(boost::accumulators::variance(accumulator) / (boost::accumulators::count(accumulator) - 1));
	double lowerBound = avg - 1.96 * stdev;
	double upperBound = avg + 1.96 * stdev;

	INFO_MESSAGE("Monte-Carlo estimate: " << avg)
	INFO_MESSAGE("Stdev: " << stdev)
	INFO_MESSAGE("95% confidence interval: [" << lowerBound << ", " << upperBound << "]")

	// calculate the theoretical value according to Black formula
	boost::math::normal s;
	RealVector sigmaVec(dim);
	(*pvol)(0, 1, sigmaVec);
	double sigma2 = blitz::sum(sigmaVec * sigmaVec);
	DEBUG_MESSAGE("sigma2: " << sigma2)
	double d1 = (log(initial(1) / K) + sigma2 * timeline->delta(0) / 2) / sqrt(sigma2 * timeline->delta(0));
	double d2 = (log(initial(1) / K) - sigma2 * timeline->delta(0) / 2) / sqrt(sigma2 * timeline->delta(0));
	double x = (initial(1) * boost::math::cdf(s, d1) - K * boost::math::cdf(s, d2)) * timeline->delta(1) / ((1 + timeline->delta(0) * initial(0)) * (1 + timeline->delta(1) * initial(1)));
	INFO_MESSAGE("theoretical: " << x)

	// check Monte-Carlo is within confidence bounds, it should be there with 95% probability
	bool withinBounds = x > lowerBound && x < upperBound;
	INFO_MESSAGE("Within confidence bounds: " << (withinBounds ? "yes" : "no"))
} // tarnpricing::caplet()

} // namespace tarnpricing
