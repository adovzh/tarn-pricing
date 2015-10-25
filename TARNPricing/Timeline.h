#ifndef __TARNPRICING_TIMELINE_H
#define __TARNPRICING_TIMELINE_H

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <blitz/array.h>
#include "Types.h"

namespace tarnpricing {

class Timeline {
public:
	// constructor is inlined
	Timeline(double start, double maturity, int periods);
	Timeline(const RealVector& points);

	// inline
	// delta(i) = T(i + 1) - T(i)
	double delta(int i) const { return deltas(i); }
	int length() const { return deltas.extent(blitz::firstDim); }
	double operator()(int i) const { return time_points(i); }

	typedef boost::shared_ptr<Timeline> Ptr;
	typedef boost::shared_ptr<const Timeline> ConstPtr;
	friend std::ostream& operator<<(std::ostream& out, const Timeline& timeline);
private:
	RealVector time_points;
	RealVector deltas;
}; // class Timeline

} // namespace tarnpricing

#endif // __TARNPRICING_TIMELINE_H