#ifndef __PARAMETERISED_VOLATILITY_H
#define __PARAMETERISED_VOLATILITY_H

#include <boost/shared_ptr.hpp>
#include "Types.h"
#include "Timeline.h"

namespace tarnpricing {

class ParameterisedVolatility {
public:
	ParameterisedVolatility(int _dim, const Timeline::Ptr _timeline, RealVector& _a, const RealVector& _c): dim(_dim), timeline(_timeline), a(_a), c(_c) {}
	int dimension() const { return dim; };
	void operator()(int i, int j, RealVector& sigma) const;

	typedef boost::shared_ptr<ParameterisedVolatility> Ptr;
	typedef boost::shared_ptr<const ParameterisedVolatility> ConstPtr;
private:
	int dim;
	Timeline::Ptr timeline;
	RealVector a;
	RealVector c;
};

inline void ParameterisedVolatility::operator()(int i, int j, RealVector& sigma) const
{
	Timeline& t = *timeline;

	sigma = c * exp(-a * (t(j) - t(i)));
}

} // namespace tarnpricing

#endif