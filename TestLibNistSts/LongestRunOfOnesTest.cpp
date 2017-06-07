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
	TEST_CLASS(LongestRunOfOnesTest)
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

		// compare each version of cusum test
		TEST_METHOD(LongestRunOfOnesCompareVersions)
		{
			auto seq = GetSeqLcg(4096);
			seq[0] = false;
			seq[1] = false;
			seq[2] = false;
			Nist::Parameters parameters;
			Nist::Results results1, results2, results3;
			Nist::Test test1(&seq, &parameters, &results1);
			Nist::Test test2(&seq, &parameters, &results2);
			Nist::Test test3(&seq, &parameters, &results3);
			LongestRunOfOnes(test1);
			LongestRunOfOnes2(test2);
			LongestRunOfOnes3(test3);
			auto& result1 = test1.GetResults().longestrunofones;
			auto& result2 = test2.GetResults().longestrunofones;
			auto& result3 = test3.GetResults().longestrunofones;
			Assert::AreEqual(result1.p_value, result2.p_value, 0.01);
			Assert::AreEqual(result1.p_value, result3.p_value, 0.01);
			Assert::AreEqual(result1.chi2, result2.chi2, 0.01);
			Assert::AreEqual(result1.chi2, result3.chi2, 0.01);
			Assert::AreEqual(result1.M, result2.M);
			Assert::AreEqual(result1.M, result3.M);
			Assert::AreEqual(result1.N, result2.N);
			Assert::AreEqual(result1.N, result3.N);
			for (int i = 0; i < result1.nu.size(); i++) {
				Assert::AreEqual(result1.nu[i], result2.nu[i]);
				Assert::AreEqual(result1.nu[i], result3.nu[i]);
			}
		}

	};
}