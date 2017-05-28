#include "stdafx.h"
#include "CppUnitTest.h"
#include <vector>
#include "../LibNistSts/BitSequence.h"
#include "../LibNistSts/Parameters.h"
#include "../LibNistSts/Test.h"
#include "../LibNistSts/common/stat_fncs.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace libNISTtest
{
	TEST_CLASS(BlockFrequencyUnitTest)
	{
	public:

		Nist::BitSequence GetSeqLcg(size_t n)
		{
			std::vector<bool> data(n);
			// make sequence of 'random'
			int lcg = 112437;
			data.resize(256);
			for (auto& bit : data)
			{
				bit = lcg & 8 == 8;
				lcg = lcg * 7612341 + 41235467;
			}
			return Nist::BitSequence(data);
		}

		// compare each version of frequency test
		TEST_METHOD(BlockFrequencyCompareVersions)
		{
			auto seq = GetSeqLcg(4096);
			seq[0] = false;
			seq[1] = false;
			seq[2] = false;
			Nist::Parameters parameters;
			Nist::Results results1, results2, results3, results4;
			Nist::Test test1(&seq, &parameters, &results1);
			Nist::Test test2(&seq, &parameters, &results2);
			Nist::Test test3(&seq, &parameters, &results3);
			Nist::Test test4(&seq, &parameters, &results4);
			BlockFrequency(test1);
			BlockFrequency2(test2);
			BlockFrequency3(test3);
			BlockFrequency4(test4);
			auto& fr1Result = test1.GetResults().blockfrequency;
			auto& fr2Result = test2.GetResults().blockfrequency;
			auto& fr3Result = test3.GetResults().blockfrequency;
			auto& fr4Result = test4.GetResults().blockfrequency;
			Assert::AreEqual(fr1Result.chi_squared, fr2Result.chi_squared, 0.01, L"Sum should be the same for v1 and v2.");
			Assert::AreEqual(fr1Result.chi_squared, fr3Result.chi_squared, 0.01, L"Sum should be the same for v1 and v3.");
			Assert::AreEqual(fr1Result.chi_squared, fr4Result.chi_squared, 0.01, L"Sum should be the same for v1 and v4.");
			Assert::AreEqual(fr1Result.p_value, fr2Result.p_value, 0.01, L"P value should be the same for v1 and v2.");
			Assert::AreEqual(fr1Result.p_value, fr3Result.p_value, 0.01, L"P value should be the same for v1 and v3.");
			Assert::AreEqual(fr1Result.p_value, fr4Result.p_value, 0.01, L"P value should be the same for v1 and v4.");
			
		}

	};
}