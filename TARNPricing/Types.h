#ifndef __TARNPRICING_TYPES_H
#define __TARNPRICING_TYPES_H

#include <boost/function.hpp>
#include <blitz/array.h>

namespace tarnpricing {

template<class T>
struct Matrix
{
	typedef blitz::Array<T,2> type;
	typedef blitz::Array<T,2>& Ref;
};

template<typename ARG>
struct BoostRNG
{
	typedef boost::function<ARG ()> type;
};

} // namespace tarnpricing

#endif // __TARNPRICING_TYPES_H