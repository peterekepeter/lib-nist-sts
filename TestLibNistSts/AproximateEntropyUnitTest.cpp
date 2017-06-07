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
	TEST_CLASS(AproximateEntropyUnitTest)
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

		// compare each version
		TEST_METHOD(AproximateEntropyCompareVersions)
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
			ApproximateEntropy(test1);
			ApproximateEntropy2(test2);
			ApproximateEntropy4(test3);
			auto& result1 = test1.GetResults().approximate_entropy;
			auto& result2 = test2.GetResults().approximate_entropy;
			auto& result3 = test3.GetResults().approximate_entropy;
			for (size_t i = 0; i < 2; i++) {
				Assert::AreEqual(result1.ApEn[i], result2.ApEn[i], 0.01);
				Assert::AreEqual(result1.ApEn[i], result3.ApEn[i], 0.01);
			}
			Assert::AreEqual(result1.P.size(), result2.P.size());
			Assert::AreEqual(result1.P.size(), result3.P.size());
			for (size_t i = 0; i < result1.P.size(); i++) {
				Assert::AreEqual(result1.P[i], result2.P[i], 0.01);
				Assert::AreEqual(result1.P[i], result3.P[i], 0.01);
			}
			Assert::AreEqual(result1.chi_squared, result2.chi_squared, 0.01);
			Assert::AreEqual(result1.chi_squared, result3.chi_squared, 0.01);
			Assert::AreEqual(result1.pp, result2.pp, 0.01);
			Assert::AreEqual(result1.pp, result3.pp, 0.01);
			Assert::AreEqual(result1.p_value, result2.p_value, 0.01);
			Assert::AreEqual(result1.p_value, result3.p_value, 0.01);
		}

	};
}