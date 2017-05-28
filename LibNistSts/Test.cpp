#include "Test.h"
#include "common/stat_fncs.h"


const Nist::Results::Frequency& Nist::Test::RunFrequency()
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

const Nist::Results::BlockFrequency& Nist::Test::RunBlockFrequency()
{
	if (this->parameters->fast == 0)
	{
		BlockFrequency(*this);
	}
	else
	{
		BlockFrequency4(*this);
	}
	return results->blockfrequency;
}

const Nist::Results& Nist::Test::RunAll()
{
	this->RunFrequency();
	this->RunBlockFrequency();
	return *results;
}
