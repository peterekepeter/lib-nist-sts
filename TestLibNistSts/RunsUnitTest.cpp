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
	TEST_CLASS(RunsUnitTest)
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
		TEST_METHOD(RunsCompareVersions)
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
			Runs(test1);
			Runs2(test2);
			Runs3(test3);
			Runs4(test4);
			auto& result1 = test1.GetResults().runs;
			auto& result2 = test2.GetResults().runs;
			auto& result3 = test3.GetResults().runs;
			auto& result4 = test4.GetResults().runs;
			Assert::AreEqual(result1.pi, result2.pi, 0.01);
			Assert::AreEqual(result1.pi, result3.pi, 0.01);
			Assert::AreEqual(result1.pi, result4.pi, 0.01);
			Assert::AreEqual(result1.p_value, result2.p_value, 0.01);
			Assert::AreEqual(result1.p_value, result3.p_value, 0.01);
			Assert::AreEqual(result1.p_value, result4.p_value, 0.01);
			Assert::AreEqual(result1.V, result2.V, 0.01);
			Assert::AreEqual(result1.V, result3.V, 0.01);
			Assert::AreEqual(result1.V, result4.V, 0.01);
			Assert::AreEqual(result1.erfc_arg, result2.erfc_arg, 0.01);
			Assert::AreEqual(result1.erfc_arg, result3.erfc_arg, 0.01);
			Assert::AreEqual(result1.erfc_arg, result4.erfc_arg, 0.01);
		}

	};
}