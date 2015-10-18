#include <iostream>
#include "Timeline.h"

namespace tarnpricing {

Timeline::Timeline(double start, double maturity, int periods): time_points(periods + 1), deltas(periods)
{
	blitz::firstIndex i;
	double dt = maturity / periods;

	time_points = start + i * dt;
	deltas = dt;
}

Timeline::Timeline(const RealVector& points): time_points(points.extent(blitz::firstDim)), deltas(points.extent(blitz::firstDim) - 1)
{
	time_points = points.copy();

	for (int i = 0; i < points.extent(blitz::firstDim) - 1; i++) deltas(i) = points(i + 1) - points(i);
}


std::ostream& operator<<(std::ostream& out, const tarnpricing::Timeline& timeline)
{
	out << timeline.time_points << std::endl;
	out << timeline.deltas << std::endl;

	return out;
}

} // namespace tarnpricing
