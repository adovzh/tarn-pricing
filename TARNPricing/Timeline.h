#ifndef __TARNPRICING_TIMELINE_H
#define __TARNPRICING_TIMELINE_H

#include <blitz/array.h>
#include <iostream>

namespace tarnpricing {

	class Timeline {
	private:
		blitz::Array<double,1> time_points;
		blitz::Array<double,1> deltas;

	public:
		// constructor is inlined
		Timeline(double start, double maturity, int periods): time_points(periods + 1), deltas(periods)
		{
			blitz::firstIndex i;
			double dt = maturity / periods;

			time_points = start + i * dt;
			deltas = dt;
		}

		// inline
		// delta(i) = T(i + 1) - T(i)
		double delta(int i) const { return deltas(i); }

		friend std::ostream& operator<<(std::ostream& out, const Timeline& timeline);
	};
}

#endif // __TARNPRICING_TIMELINE_H