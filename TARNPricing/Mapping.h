#ifndef __TARNPRICING_MAPPING_H
#define __TARNPRICING_MAPPING_H

template<typename price_process, typename mapped_rv>
class Mapping
{
private:
	price_process& process;

public:
	Mapping();
	double mapping(mapped_rv x);
};

template<typename price_process, typename mapped_rv>
Mapping<price_process, mapped_rv>::Mapping(price_process& _process): process(_process)
{
	// set timeline from payoff to process
}

template<typename price_process, typename mapped_rv>
Mapping<price_process, mapped_rv>::mapping(mapped_rv x)
{
	// implement the mapping
	return 0;
}

#endif // __TARNPRICING_MAPPING_H