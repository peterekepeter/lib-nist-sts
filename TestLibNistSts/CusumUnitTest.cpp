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
	TEST_CLASS(CumulativeSumsTest)
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
		TEST_METHOD(CusumCompareVersions)
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
			CumulativeSums(test1);
			CumulativeSums2(test2);
			CumulativeSums3(test3);
			auto& result1 = test1.GetResults().cusum;
			auto& result2 = test2.GetResults().cusum;
			auto& result3 = test3.GetResults().cusum;
			Assert::AreEqual(result1.p_valueA, result2.p_valueA, 0.01);
			Assert::AreEqual(result1.p_valueA, result3.p_valueA, 0.01);
			Assert::AreEqual(result1.p_valueB, result2.p_valueB, 0.01);
			Assert::AreEqual(result1.p_valueB, result3.p_valueB, 0.01);
			Assert::AreEqual(result1.sum1A, result2.sum1A, 0.01);
			Assert::AreEqual(result1.sum2A, result2.sum2A, 0.01);
			Assert::AreEqual(result1.sum1B, result2.sum1B, 0.01);
			Assert::AreEqual(result1.sum2B, result2.sum2B, 0.01);
			Assert::AreEqual(result1.z, result2.z, 0.01);
			Assert::AreEqual(result1.z, result2.z, 0.01);
			Assert::AreEqual(result1.zrev, result2.zrev, 0.01);
			Assert::AreEqual(result1.zrev, result2.zrev, 0.01);
		}

	};
}