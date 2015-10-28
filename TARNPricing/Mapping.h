#ifndef __TARNPRICING_MAPPING_H
#define __TARNPRICING_MAPPING_H

#include "Types.h"
#include "Payoff.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

/**
 * Template that maps a set of random numbers to a stochastic model, and then 
 * to a single-dimensional payoff.
 *
 * Template parameter MODEL must implement:
 *   - member function timeline(), which sets a timeline;
 *   - member function initial(), which returns an initial value of an underlying
 *   - operator()(mapped_rv, LowerTriangularMatrix&), which fills the specified matrix with the simulated values
 */
template<typename MODEL, typename mapped_rv, typename MODEL_PTR = MODEL::Ptr>
class Mapping
{
public:
	/// Constructor
	Mapping(const MODEL_PTR& _model, const Payoff::ConstPtr& _payoff);
	/// Mapping of random numbers to a discounted payoff
	double mapping(mapped_rv x);
private:
	const MODEL_PTR model; ///< stochastic model
	const Payoff::ConstPtr payoff; ///< smart pointer to a payoff
	LowerTriangularMatrix rates; ///< auxiliary matrix to pass simulated rates from the model to the payoff functor
};

template<typename MODEL, typename mapped_rv, typename MODEL_PTR>
Mapping<MODEL, mapped_rv, MODEL_PTR>::Mapping(const MODEL_PTR& _model, const Payoff::ConstPtr& _payoff): model(_model), payoff(_payoff), rates(_payoff->timeline()->length())
{
	// set timeline from payoff to model
	model->timeline(payoff->timeline());
	rates.setColumn(0, model->initial());
	// rates(blitz::Range::all(), 0) = model->initial();
}

template<typename MODEL, typename mapped_rv, typename MODEL_PTR>
double Mapping<MODEL, mapped_rv, MODEL_PTR>::mapping(mapped_rv x)
{
	// implement the mapping
	(*model)(x, rates);
	return (*payoff)(rates);
}

} // namespace tarnpricing

#endif // __TARNPRICING_MAPPING_H