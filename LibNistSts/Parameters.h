#pragma once
#include <stddef.h>
#include <limits>

namespace Nist
{
	class Parameters
	{
	public:
		size_t  n = std::numeric_limits<size_t>::max(); // default: = std::numeric_limits<size_t>.max()
		int		blockFrequencyBlockLength = 128;
		int		nonOverlappingTemplateBlockLength = 9;
		int		overlappingTemplateBlockLength = 9;
		int		serialBlockLength = 16;
		int		linearComplexitySequenceLength = 500;
		int		approximateEntropyBlockLength = 10;
		int		numOfBitStreams = 1;
		// New Flag
		int		fast = 1;
	};
}
