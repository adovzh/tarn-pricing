#ifndef __PARAMETERISED_VOLATILITY_H
#define __PARAMETERISED_VOLATILITY_H

#include <boost/shared_ptr.hpp>
#include "Types.h"
#include "Timeline.h"

namespace tarnpricing {

/**
 * Covariance matrix structure at T_i for L_j, specified as sigma(i,j) = c(i,j) * exp(-a(i,j) * abs(T_j - T_i)).
 */
class ParameterisedVolatility {
public:
	/// Constructor
	ParameterisedVolatility(int _dim, const Timeline::ConstPtr& _timeline, RealVector& _a, const RealVector& _c): dim(_dim), timeline(_timeline), a(_a), c(_c) {}
	/// Dimension of sigma(i,j), a and c vectors.
	int dimension() const { return dim; };
	/// Put sigma(i, j) toa specified vector sigma.
	void operator()(int i, int j, RealVector& sigma) const;

	typedef boost::shared_ptr<ParameterisedVolatility> Ptr;
	typedef boost::shared_ptr<const ParameterisedVolatility> ConstPtr;
private:
	int dim; ///< dimension of sigma(i,j), a and c vectors.
	const Timeline::ConstPtr timeline; ///< timeline provided
	RealVector a; ///< parameters, see class doc comment
	RealVector c; ///< parameters, see class doc comment
};

inline void ParameterisedVolatility::operator()(int i, int j, RealVector& sigma) const
{
	const Timeline& t = *timeline;

	sigma = c * exp(-a * (t(j) - t(i)));
}

} // namespace tarnpricing

#endif