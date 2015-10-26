#include <algorithm>
#include <boost/make_shared.hpp>
#include "CapletPayoff.h"
#include "Logging.h"

namespace tarnpricing {

CapletPayoff::CapletPayoff(double start, double fixing, double payment, double strike): K(strike)
{
	RealVector vTimeline(3);
	vTimeline = start, fixing, payment;
	Timeline::ConstPtr timeline = boost::make_shared<const Timeline>(vTimeline);
	setTimeline(timeline);
}

double CapletPayoff::operator()(const LowerTriangularMatrix& underlying) const
{
	double payoff = std::max((underlying(1, 1) - K) * timeline()->delta(1), 0.0);
	double numeraire = (1 + underlying(0, 0) * timeline()->delta(0)) * (1 + underlying(1, 1) * timeline()->delta(1));

	TRACE_MESSAGE("L(1,1) = " << underlying(1, 1) << " against strike = " << K)
	TRACE_MESSAGE("discounted payoff is " << payoff / numeraire)

	return payoff / numeraire;
}

} // namespace tarnpricing