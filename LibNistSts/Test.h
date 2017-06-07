#pragma once
#include "BitSequence.h"
#include "Parameters.h"
#include "Results.h"

namespace Nist
{
	// represents a test session, includes input, parameters and output
	class Test
	{
	public:
		Test(BitSequence* bit_stream, Parameters* parameters = nullptr, Results* results = nullptr)
			: bitStream(bit_stream),
			  parameters(parameters),
			  results(results),
			  ownedParameters(false),
			  ownedResults(false)
		{
			if (parameters == nullptr)
			{
				// use defaults
				this->parameters = new Parameters();
				ownedParameters = true;
			}
			if (results == nullptr)
			{
				// create struct
				this->results = new Results();
				ownedResults = true;
			}
		}

		~Test()
		{
			if (ownedParameters && parameters != nullptr)
			{
				delete parameters;
				ownedParameters = false;
				parameters = nullptr;
			}
			if (ownedResults && results != nullptr)
			{
				delete results;
				ownedResults = false;
				results = nullptr;
			}
		}


		// get sparse array of bits, each byte is one bit, used by original NIST implementation
		const unsigned char* GetEpsilon() const
		{
			return bitStream->GetEpsilon();
		}

		// get packed array of bits, all 8 bits of each byte is used, used by fast NIST implementation
		const unsigned char* GetArray() const
		{
			return bitStream->GetArray();
		}

		// get the number of bits, in case of test it's the minimum between bit stream length and parameters->n
		size_t GetN() const
		{
			size_t n1 = parameters->n;
			size_t n2 = bitStream->GetN();
			return n1 < n2 ? n1 : n2;
		}

		const Results& GetResults() const
		{
			return *results;
		}

		Results& GetWriteableResults()
		{
			return *results;
		}

		const Parameters& GetParameters() const
		{
			return *parameters;
		}

		const Results::Frequency& RunTestFrequency();

		const Results::BlockFrequency& RunTestBlockFrequency();

		const Results::Cusum& RunTestCumulativeSums();

		const Results::Runs& RunTestRuns();

		const Results::LongestRunOfOnes& RunTestLongestRunOfOnes();

		const Results::Rank& RunTestRank();

		const Results::Dft& RunTestDiscreteFourierTransform();

		const Results::NonOverlapping& RunTestNonOverlappingTemplateMatchings();

		const Results::Overlapping& RunTestOverlappingTemplateMatchings();

		const Results::Universal& RunTestUniversal();

		const Results::ApproximateEntropy& RunTestApproximateEntropy();

		const Results::RandomExcursion& RunTestRandomExcursions();

		const Results::RandomExcursionVariant& RunTestRandomExcursionsVariant();

		const Results::LinearComplexity& RunTestLinearComplexity();

		const Results::Serial& RunTestSerial();

		const Results& RunTestAll();

	private:
		BitSequence* bitStream;
		Parameters* parameters;
		Results* results;
		bool ownedParameters;
		bool ownedResults;
	};
}
