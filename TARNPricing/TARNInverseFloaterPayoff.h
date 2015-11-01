#ifndef __TARN_INVERSE_FLOATER_PAYOFF_H
#define __TARN_INVERSE_FLOATER_PAYOFF_H

#include <boost/function.hpp>

#include "Payoff.h"

namespace tarnpricing {

class TARNInverseFloaterPayoff: public Payoff
{
public:
	TARNInverseFloaterPayoff(const Timeline::ConstPtr& _timeline, double _targetCoupon, double _c1, double _K, double _alpha);
	double operator()(const LowerTriangularMatrix& underlying) const;
private:
	double targetCoupon;
	double c1;
	double K;
	double alpha;
};

}

#endif