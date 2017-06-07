#include "Test.h"
#include "common/stat_fncs.h"


const Nist::Results::Frequency& Nist::Test::RunTestFrequency()
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

const Nist::Results::BlockFrequency& Nist::Test::RunTestBlockFrequency()
{
	if (this->parameters->fast == 0)
	{
		BlockFrequency_v1(*this);
	}
	else
	{
		BlockFrequency_v2(*this);
	}
	return results->blockfrequency;
}

const Nist::Results::Cusum & Nist::Test::RunTestCumulativeSums()
{
	if (this->parameters->fast == 0)
	{
		CumulativeSums_v1(*this);
	}
	else
	{
		CumulativeSums_v2(*this);
	}
	return results->cusum;
}

const Nist::Results::Runs & Nist::Test::RunTestRuns()
{
	if (this->parameters->fast == 0)
	{
		Runs_v1(*this);
	}
	else
	{
		Runs_v2(*this);
	}
	return results->runs;
}

const Nist::Results::LongestRunOfOnes & Nist::Test::RunTestLongestRunOfOnes()
{
	if (this->parameters->fast == 0)
	{
		LongestRunOfOnes_v1(*this);
	}
	else
	{
		LongestRunOfOnes_v2(*this);
	}
	return results->longestrunofones;
}

const Nist::Results::Rank & Nist::Test::RunTestRank()
{
	if (this->parameters->fast == 0)
	{
		Rank_v1(*this);
	}
	else
	{
		Rank_v2(*this);
	}
	return results->rank;
}

const Nist::Results::Dft & Nist::Test::RunTestDiscreteFourierTransform()
{
	if (this->parameters->fast == 0)
	{
		DiscreteFourierTransform_v1(*this);
	}
	else
	{
		DiscreteFourierTransform_v2(*this);
	}
	return results->dft;
}

const Nist::Results::NonOverlapping & Nist::Test::RunTestNonOverlappingTemplateMatchings()
{
	if (this->parameters->fast == 0)
	{
		NonOverlappingTemplateMatchings_v1(*this);
	}
	else
	{
		NonOverlappingTemplateMatchings_v2(*this);
	}
	return results->nonoverlapping;
}

const Nist::Results::Overlapping & Nist::Test::RunTestOverlappingTemplateMatchings()
{
	if (this->parameters->fast == 0)
	{
		OverlappingTemplateMatchings_v1(*this);
	}
	else
	{
		OverlappingTemplateMatchings_v2(*this);
	}
	return results->overlapping;
}

const Nist::Results::Universal & Nist::Test::RunTestUniversal()
{
	if (this->parameters->fast == 0)
	{
		Universal_v1(*this);
	}
	else
	{
		Universal_v2(*this);
	}
	return results->universal;
}

const Nist::Results::ApproximateEntropy & Nist::Test::RunTestApproximateEntropy()
{
	if (this->parameters->fast == 0)
	{
		ApproximateEntropy_v1(*this);
	}
	else
	{
		ApproximateEntropy_v2(*this);
	}
	return results->approximate_entropy;
}

const Nist::Results::RandomExcursion & Nist::Test::RunTestRandomExcursions()
{
	if (this->parameters->fast == 0)
	{
		RandomExcursions_v1(*this);
	}
	else
	{
		RandomExcursions_v2(*this);
	}
	return results->random_excursion;
}

const Nist::Results::RandomExcursionVariant & Nist::Test::RunTestRandomExcursionsVariant()
{
	if (this->parameters->fast == 0)
	{
		RandomExcursionsVariant_v1(*this);
	}
	else
	{
		RandomExcursionsVariant_v2(*this);
	}
	return results->random_excursion_variant;
}

const Nist::Results::LinearComplexity & Nist::Test::RunTestLinearComplexity()
{
	if (this->parameters->fast == 0)
	{
		LinearComplexity_v1(*this);
	}
	else
	{
		LinearComplexity_v2(*this);
	}
	return results->linear_complexity;
}

const Nist::Results::Serial & Nist::Test::RunTestSerial()
{
	if (this->parameters->fast == 0)
	{
		Serial_v1(*this);
	}
	else
	{
		Serial_v2(*this);
	}
	return results->serial;
}

const Nist::Results& Nist::Test::RunTestAll()
{
	this->RunTestFrequency();
	this->RunTestBlockFrequency();
	this->RunTestCumulativeSums();
	this->RunTestRuns();
	this->RunTestLongestRunOfOnes();
	this->RunTestRank();
	this->RunTestDiscreteFourierTransform();
	this->RunTestNonOverlappingTemplateMatchings();
	this->RunTestOverlappingTemplateMatchings();
	this->RunTestUniversal();
	this->RunTestApproximateEntropy();
	this->RunTestRandomExcursions();
	this->RunTestRandomExcursionsVariant();
	this->RunTestLinearComplexity();
	this->RunTestSerial();
	return *results;
}
