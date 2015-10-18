#ifndef __TARNPRICING_TIMELINE_H
#define __TARNPRICING_TIMELINE_H

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <blitz/array.h>
#include "Types.h"

namespace tarnpricing {

	class Timeline {
	private:
		RealVector time_points;
		RealVector deltas;

		typedef boost::shared_ptr<Timeline> Ptr;
	public:
		// constructor is inlined
		Timeline(double start, double maturity, int periods);
		Timeline(const RealVector& points);

		// inline
		// delta(i) = T(i + 1) - T(i)
		double delta(int i) const { return deltas(i); }

		friend std::ostream& operator<<(std::ostream& out, const Timeline& timeline);
	};
}

#endif // __TARNPRICING_TIMELINE_H