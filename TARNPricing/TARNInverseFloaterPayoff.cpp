#include "TARNInverseFloaterPayoff.h"
#include "Logging.h"

namespace tarnpricing {

TARNInverseFloaterPayoff::TARNInverseFloaterPayoff(const Timeline::ConstPtr& _timeline, double _targetCoupon, double _c1, double _K, double _alpha):
	Payoff(_timeline), targetCoupon(_targetCoupon), c1(_c1), K(_K), alpha(_alpha)
{
}

double TARNInverseFloaterPayoff::operator()(const LowerTriangularMatrix& underlying) const
{
	const Timeline::ConstPtr& t = timeline();
	double numeraire = (1 + underlying(0, 0) * t->delta(0)) * (1 + underlying(1, 1) * t->delta(1));
	double couponSum = c1 * t->delta(1);
	double payoff = couponSum / numeraire;
	int i = 2;
	double c; // coupon
	
	TRACE_MESSAGE("T_1: " << couponSum)

	while (i < t->length())
	{
		double u = underlying(i, i);
		double d = t->delta(i);

		c = std::max(K - alpha * underlying(i, i) * t->delta(i), 0.0);
		TRACE_MESSAGE("T_" << i <<": " << c)

		numeraire *= 1 + underlying(i, i) * t->delta(i);
		couponSum += c;

		if (couponSum > targetCoupon)	
		{
			c -= couponSum - targetCoupon;
			payoff += c / numeraire;
			TRACE_MESSAGE("T_" << i <<" (actual): " << c)
			break;
		}
		else payoff += c / numeraire;

		i++;
	}

	if (couponSum < targetCoupon)
	{
		// target coupon guarantee
		c = targetCoupon - couponSum;
		payoff += c / numeraire;
		TRACE_MESSAGE("Target coupon guarantee: " << c)
	}

	// redemption
	payoff += 1 / numeraire;

	TRACE_MESSAGE("payoff: " << payoff)

	return payoff;
}

}