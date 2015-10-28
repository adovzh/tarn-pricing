#ifndef __PAYOFF_H
#define __PAYOFF_H

#include <boost/shared_ptr.hpp>
#include "Timeline.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

/**
 * Abstract class representing a payoff of an instrument. Subclasses should implement operator()
 * and provide their specific payoff implementations.
 */
class Payoff
{
public:
	Payoff(const Timeline::ConstPtr& _timeline): pTimeline(_timeline) {}
	/// obtain a timeline, used in this payoff
	const Timeline::ConstPtr& timeline() const { return pTimeline; }
	/// fure virtual function, representing payoff
	virtual double operator()(const LowerTriangularMatrix& underlying) const = 0;

	// type definitions
	typedef boost::shared_ptr<Payoff> Ptr;
	typedef boost::shared_ptr<const Payoff> ConstPtr;
protected:
	/// subsclasses are allowed to create a payoff without specifying a timeline
	Payoff() {}
	/// set timeline later
	void setTimeline(const Timeline::ConstPtr& _timeline) { pTimeline = _timeline; }
private:
	Timeline::ConstPtr pTimeline; ///> timeline, used in this payoff
};

} // namespace tarnpricing

#endif // __PAYOFF_H