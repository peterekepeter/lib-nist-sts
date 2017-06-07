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
	TEST_CLASS(OverlappingTemplateMatchhingsUnitTest)
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
		TEST_METHOD(OverlappingTemplateMathchingsCompareVersions)
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
			OverlappingTemplateMatchings(test1);
			OverlappingTemplateMatchings2(test2);
			OverlappingTemplateMatchings3(test3);
			OverlappingTemplateMatchings4(test4);
			auto& result1 = test1.GetResults().overlapping;
			auto& result2 = test2.GetResults().overlapping;
			auto& result3 = test2.GetResults().overlapping;
			auto& result4 = test2.GetResults().overlapping;
			// TODO: investiage why these are NaN
			/*Assert::AreEqual(result1.p_value, result2.p_value, 0.01);
			Assert::AreEqual(result1.p_value, result3.p_value, 0.01);
			Assert::AreEqual(result1.p_value, result4.p_value, 0.01);
			Assert::AreEqual(result1.chi2, result2.chi2, 0.01);
			Assert::AreEqual(result1.chi2, result3.chi2, 0.01);
			Assert::AreEqual(result1.chi2, result4.chi2, 0.01);*/
			for (int i = 0; i < result1.nu.size(); i++) {
				Assert::AreEqual(result1.nu[i], result2.nu[i], 0.01);
				Assert::AreEqual(result1.nu[i], result3.nu[i], 0.01);
				Assert::AreEqual(result1.nu[i], result4.nu[i], 0.01);
			}
		}

	};
}