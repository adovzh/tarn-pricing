#ifndef __TARNPRICING_TYPES_H
#define __TARNPRICING_TYPES_H

#include <boost/function.hpp>
#include <blitz/array.h>
#include "TARNPricing.h"

TP_NAMESPACE(tarnpricing)

template<typename T>
struct Vector
{
	typedef blitz::Array<T,1> type;
};

typedef blitz::Array<double,1> RealVector;

template<typename T>
struct Matrix
{
	typedef blitz::Array<T,2> type;
	typedef blitz::Array<T,2>& Ref;
};

typedef blitz::Array<double,2> RealMatrix;

template<typename ARG>
struct BoostRNG
{
	typedef boost::function<ARG ()> type;
};

typedef boost::function<RealVector (double, double)> Volatility;

TP_NAMESPACE_END

#endif // __TARNPRICING_TYPES_H