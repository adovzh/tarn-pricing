#ifndef __PAYOFF_H
#define __PAYOFF_H

#include <boost/shared_ptr.hpp>
#include "Timeline.h"
#include "LowerTriangularMatrix.h"

namespace tarnpricing {

class Payoff
{
public:
	Payoff(const Timeline::ConstPtr& _timeline): pTimeline(_timeline) {}
	const Timeline::ConstPtr& timeline() const { return pTimeline; }

	virtual double operator()(const LowerTriangularMatrix& underlying) const = 0;

	// type definitions
	typedef boost::shared_ptr<Payoff> Ptr;
	typedef boost::shared_ptr<const Payoff> ConstPtr;
protected:
	Payoff() {}
	void setTimeline(const Timeline::ConstPtr& _timeline) { pTimeline = _timeline; }
private:
	Timeline::ConstPtr pTimeline;
};

} // namespace tarnpricing

#endif // __PAYOFF_H