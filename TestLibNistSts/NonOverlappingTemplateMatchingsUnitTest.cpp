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
	TEST_CLASS(NonOverlappingTemplateMatchhingsUnitTest)
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
		TEST_METHOD(NonOverlappingTemplateMathchingsCompareVersions)
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
			NonOverlappingTemplateMatchings(test1);
			NonOverlappingTemplateMatchings2(test2);
			NonOverlappingTemplateMatchings4(test3);
			auto& result1 = test1.GetResults().nonoverlapping;
			auto& result2 = test2.GetResults().nonoverlapping;
			auto& result3 = test2.GetResults().nonoverlapping;

			Assert::AreEqual(result1.templates, result2.templates);
			Assert::AreEqual(result1.templates, result3.templates);
			Assert::AreEqual(result1.chi2.size(), result2.chi2.size());
			Assert::AreEqual(result1.chi2.size(), result3.chi2.size());
			Assert::AreEqual(result1.p_value.size(), result2.p_value.size());
			Assert::AreEqual(result1.p_value.size(), result3.p_value.size());
			Assert::AreEqual(result1.W.size(), result2.W.size());
			Assert::AreEqual(result1.W.size(), result3.W.size());

			for (size_t i = 0; i < result1.chi2.size(); i++) {
				Assert::AreEqual(result1.chi2[i], result2.chi2[i], 0.01);
				Assert::AreEqual(result1.chi2[i], result3.chi2[i], 0.01);
			}
			for (size_t i = 0; i < result1.p_value.size(); i++) {
				Assert::AreEqual(result1.p_value[i], result2.p_value[i], 0.01);
				Assert::AreEqual(result1.p_value[i], result3.p_value[i], 0.01);
			}
			for (size_t i = 0; i < result1.W.size(); i++) {
				Assert::AreEqual(result1.W[i], result2.W[i], 0.01);
				Assert::AreEqual(result1.W[i], result3.W[i], 0.01);
			}
		}

	};
}