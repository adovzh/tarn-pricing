#ifndef __TARNPRICING_MAPPING_H
#define __TARNPRICING_MAPPING_H

#include "Types.h"
#include "Payoff.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

template<typename MODEL, typename mapped_rv, typename MODEL_PTR = MODEL::Ptr>
class Mapping
{
public:
	Mapping(const MODEL_PTR& _model, const Payoff::ConstPtr& _payoff);
	double mapping(mapped_rv x);
private:
	const MODEL_PTR model;
	const Payoff::ConstPtr payoff;
	LowerTriangularMatrix rates;
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