#ifndef __CAPLET_PAYOFF_H
#define __CAPLET_PAYOFF_H

#include "Types.h"
#include "Payoff.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

class CapletPayoff: public Payoff
{
public:
	CapletPayoff(double start, double fixing, double payment, double strike);
	double operator()(const LowerTriangularMatrix& underlying) const;
private:
	const double K;
};

} // tarnpricing

#endif // __CAPLET_PAYOFF