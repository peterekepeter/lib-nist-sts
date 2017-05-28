#include "Test.h"
#include "common/stat_fncs.h"


Nist::Results::Frequency& Nist::Test::RunFrequency()
{
	if (this->parameters->fast == 0)
	{
		Frequency_v1(*this);
	} 
	else
	{
		Frequency_v2(*this);
	}
	return results->frequency;
}
