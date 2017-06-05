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
	TEST_CLASS(RankUnitTest)
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
		TEST_METHOD(RankCompareVersions)
		{
			auto seq = GetSeqLcg(4096);
			seq[0] = false;
			seq[1] = false;
			seq[2] = false;
			Nist::Parameters parameters;
			Nist::Results results1, results2;
			Nist::Test test1(&seq, &parameters, &results1);
			Nist::Test test2(&seq, &parameters, &results2);
			Rank(test1);
			Rank2(test2);
			auto& result1 = test1.GetResults().rank;
			auto& result2 = test2.GetResults().rank;
			Assert::AreEqual(result1.p_value, result2.p_value, 0.01);
			Assert::AreEqual(result1.chi_squared, result2.chi_squared, 0.01);
			Assert::AreEqual(result1.N, result2.N, 0.01);
			Assert::AreEqual(result1.F_30, result2.F_30, 0.01);
			Assert::AreEqual(result1.F_31, result2.F_31, 0.01);
			Assert::AreEqual(result1.F_32, result2.F_32, 0.01);
			Assert::AreEqual(result1.p_30, result2.p_30, 0.01);
			Assert::AreEqual(result1.p_31, result2.p_31, 0.01);
			Assert::AreEqual(result1.p_32, result2.p_32, 0.01);
		}

	};
}