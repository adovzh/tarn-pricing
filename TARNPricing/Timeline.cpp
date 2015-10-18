#include <iostream>
#include "Timeline.h"

namespace tarnpricing {

	std::ostream& operator<<(std::ostream& out, const tarnpricing::Timeline& timeline)
	{
		out << timeline.time_points << std::endl;
		out << timeline.deltas << std::endl;

		return out;
	}

}
