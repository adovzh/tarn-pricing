#ifndef __CAPLET_PAYOFF_H
#define __CAPLET_PAYOFF_H

#include "Types.h"
#include "Payoff.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

class CapletPayoff: public Payoff
{
public:
	/// Construct a payoff of a caplet from a start date, fixing date and a payment date
	CapletPayoff(double start, double fixing, double payment, double strike);
	/// Calculate the caplet payoff, given the realisation of the underlying
	double operator()(const LowerTriangularMatrix& underlying) const;
private:
	const double K; ///< Caplet strike
};

} // tarnpricing

#endif // __CAPLET_PAYOFF