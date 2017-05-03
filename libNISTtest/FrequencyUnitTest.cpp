#include "stdafx.h"
#include "CppUnitTest.h"
#include <vector>
#include "../libNIST/BitSequence.h"
#include "../libNIST/Parameters.h"
#include "../libNIST/Test.h"
#include "../libNIST/stat_fncs.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace libNISTtest
{
	TEST_CLASS(FrequencyUnitTest)
	{
	public:

		Nist::BitSequence GetSeqOsc(size_t n)
		{
			std::vector<bool> data(n);
			bool osc = false;
			for (auto& bit : data)
			{
				bit = osc;
				osc = !osc;
			}
			return Nist::BitSequence(data);
		}

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

		TEST_METHOD(FrequencySlowSeq8)
		{
			auto seq8 = GetSeqOsc(8);
			Nist::Parameters paramsA;
			
			paramsA.fast = 0;
			Nist::Test test(&seq8, &paramsA);

			auto results = test.RunFrequency();
			Assert::IsTrue(results.sum == 0.0, L"Sum should be 0");
		}

		TEST_METHOD(FrequencyFastSeq8)
		{
			auto seq8 = GetSeqOsc(8);
			std::vector<bool> data(8);

			Nist::Parameters paramsA;

			paramsA.fast = 1;
			Nist::Test test(&seq8, &paramsA);

			auto results = test.RunFrequency();
			Assert::IsTrue(results.sum == 0.0, L"Sum should be 0");
		}


		TEST_METHOD(FrequencySlowSeq256)
		{
			auto seq = GetSeqOsc(256);
			Nist::Parameters paramsA;

			paramsA.fast = 0;
			Nist::Test test(&seq, &paramsA);

			auto results = test.RunFrequency();
			Assert::IsTrue(results.sum == 0.0, L"Sum should be 0");
		}

		TEST_METHOD(FrequencyFastSeq256)
		{
			auto seq = GetSeqOsc(256);
			Nist::Parameters paramsA;

			paramsA.fast = 1;
			Nist::Test test(&seq, &paramsA);

			auto results = test.RunFrequency();
			Assert::IsTrue(results.sum == 0.0, L"Sum should be 0");
		}

		// compare each version of frequency test
		TEST_METHOD(FrequencyCompareVersions)
		{
			auto seq = GetSeqLcg(256);
			Nist::Parameters parameters;
			Nist::Results results1, results2, results3, results4;
			Nist::Test test1(&seq, &parameters, &results1);
			Nist::Test test2(&seq, &parameters, &results2);
			Nist::Test test3(&seq, &parameters, &results3);
			Nist::Test test4(&seq, &parameters, &results4);
			Frequency(test1);
			Frequency2(test2);
			Frequency3(test3);
			Frequency4(test4);
			auto& fr1Result = test1.GetResults().frequency;
			auto& fr2Result = test2.GetResults().frequency;
			auto& fr3Result = test3.GetResults().frequency;
			auto& fr4Result = test4.GetResults().frequency;
			Assert::AreEqual(fr1Result.sum, fr2Result.sum, 0.01, L"Sum should be the same for v1 and v2.");
			Assert::AreEqual(fr1Result.sum, fr3Result.sum, 0.01, L"Sum should be the same for v1 and v3.");
			Assert::AreEqual(fr1Result.sum, fr4Result.sum, 0.01, L"Sum should be the same for v1 and v4.");
			Assert::AreEqual(fr1Result.p_value, fr2Result.p_value, 0.01, L"P value should be the same for v1 and v2.");
			Assert::AreEqual(fr1Result.p_value, fr3Result.p_value, 0.01, L"P value should be the same for v1 and v3.");
			Assert::AreEqual(fr1Result.p_value, fr4Result.p_value, 0.01, L"P value should be the same for v1 and v4.");
			Assert::AreEqual(fr1Result.sum_n, fr2Result.sum_n, 0.01, L"Sum N value should be the same for v1 and v2.");
			Assert::AreEqual(fr1Result.sum_n, fr2Result.sum_n, 0.01, L"Sum N value should be the same for v1 and v3.");
			Assert::AreEqual(fr1Result.sum_n, fr2Result.sum_n, 0.01, L"Sum N value should be the same for v1 and v4.");
		}

	};
}